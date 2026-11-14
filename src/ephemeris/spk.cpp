/**********************************************************************
 * spk.cpp
 *
 * Copyright 2024 W. Greg Visotski
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 *         limitations under the License.
 **********************************************************************/

#include "astro/ephemeris/spk.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <mutex>

#include "astro/core/error.h"

namespace astro {

namespace {

constexpr std::size_t kRecordBytes = 1024;
constexpr std::size_t kMaxChain    = 32; // longest center chain followed (real kernels need 2-3 hops)
constexpr int         kFrameJ2000  = 1;  // NAIF frame id "J2000": the ICRF axes in the JPL planetary kernels

template <class T>
T swapBytes(T v) {
    unsigned char b[sizeof(T)];
    std::memcpy(b, &v, sizeof(T));
    std::reverse(b, b + sizeof(T));
    std::memcpy(&v, b, sizeof(T));
    return v;
}

struct Segment {
    SpkSegmentInfo info;
    std::int64_t   startAddr{0}; // 1-based double-word address
    std::int64_t   endAddr{0};
    double         init{0.0};
    double         intlen{0.0};
    std::int64_t   rsize{0};
    std::int64_t   nrec{0};
    // Cache of the most recently used record.
    mutable std::int64_t        cachedIndex{-1};
    mutable std::vector<double> cache;
};

} // namespace

struct SpkEphemeris::Impl {
    std::string                 path;
    mutable std::ifstream       file;
    mutable std::mutex          mutex;
    bool                        swap{false};
    bool                        big{false};
    std::vector<Segment>        segs;
    std::vector<SpkSegmentInfo> infos;

    template <class T>
    T fix(T v) const {
        return swap ? swapBytes(v) : v;
    }

    void readBytes(std::int64_t offset, void* dst, std::size_t n) const {
        file.clear();
        file.seekg(offset);
        file.read(static_cast<char*>(dst), static_cast<std::streamsize>(n));
        if (!file) {
            throw Error("SPK: read failed in " + path);
        }
    }

    void readDoubles(std::int64_t addr, std::size_t n, double* out) const {
        readBytes((addr - 1) * 8, out, n * sizeof(double));
        if (swap) {
            for (std::size_t i = 0; i < n; ++i) {
                out[i] = swapBytes(out[i]);
            }
        }
    }

    /// Read and check the type 2/3 directory (INIT, INTLEN, RSIZE, N) at the end of a segment:
    /// the N records of RSIZE doubles must exactly fill the segment and lie inside the file.
    void readDirectory(Segment& seg, std::int64_t fileDoubles) const {
        const std::string where = " in segment '" + seg.info.name + "' of " + path;
        if (seg.startAddr < 1 || seg.endAddr < seg.startAddr + 4 || seg.endAddr > fileDoubles) {
            throw Error("SPK: segment addresses out of range" + where);
        }
        double dir[4];
        readDoubles(seg.endAddr - 3, 4, dir);
        const double span  = static_cast<double>(seg.endAddr - seg.startAddr - 3); // doubles before dir
        const int    ncomp = seg.info.type == 2 ? 3 : 6;
        const auto   isInt = [](double x) { return std::isfinite(x) && x == std::floor(x); };
        if (!(std::isfinite(dir[0]) && std::isfinite(dir[1]) && dir[1] > 0.0)
            || !(isInt(dir[2]) && dir[2] >= 2.0 + ncomp && dir[2] <= span)
            || !(isInt(dir[3]) && dir[3] >= 1.0 && dir[2] * dir[3] == span) || std::fmod(dir[2] - 2.0, ncomp) != 0.0) {
            throw Error("SPK: corrupt segment directory" + where);
        }
        seg.init   = dir[0];
        seg.intlen = dir[1];
        seg.rsize  = static_cast<std::int64_t>(dir[2]);
        seg.nrec   = static_cast<std::int64_t>(dir[3]);
    }

    void open() {
        file.open(path, std::ios::binary);
        if (!file) {
            throw Error("SPK: cannot open " + path);
        }
        unsigned char rec[kRecordBytes];
        readBytes(0, rec, kRecordBytes);
        // The id word is "DAF/SPK ", or "NAIF/DAF" in kernels written by early SPICE Toolkits (the
        // de405.bsp that NAIF distributes is one).
        const bool older = std::memcmp(rec, "NAIF/DAF", 8) == 0;
        if (!older && std::memcmp(rec, "DAF/SPK ", 8) != 0) {
            throw Error("SPK: not a DAF/SPK file: " + path);
        }
        std::int32_t nd, ni, fward;
        std::memcpy(&nd, rec + 8, 4);
        std::memcpy(&ni, rec + 12, 4);
        std::memcpy(&fward, rec + 76, 4);
        const bool        hostBig = std::endian::native == std::endian::big;
        const std::string fmt(reinterpret_cast<const char*>(rec) + 88, 8);
        if (fmt == "BIG-IEEE") {
            big = true;
        } else if (fmt == "LTL-IEEE") {
            big = false;
        } else if (older) {
            // Those kernels were also written before the file record had a binary-format word (their
            // bytes 88-95 are zero). An SPK file has ND = 2 and NI = 6, which read as such in one byte
            // order only, so that is the order of the file; jplephem finds it the same way. If
            // neither order fits, the ND/NI check below rejects the file.
            big = (nd == 2 && ni == 6) == hostBig;
        } else {
            throw Error("SPK: unsupported binary format '" + fmt + "'");
        }
        swap  = big != hostBig;
        nd    = fix(nd);
        ni    = fix(ni);
        fward = fix(fward);
        if (nd != 2 || ni != 6) {
            throw Error("SPK: unexpected ND/NI in " + path);
        }
        const int ss = nd + (ni + 1) / 2; // summary size in doubles
        // A summary record holds 3 control words and then the summaries.
        const int maxSummaries = static_cast<int>(kRecordBytes / 8 - 3) / ss;

        file.seekg(0, std::ios::end);
        const std::int64_t fileRecords = static_cast<std::int64_t>(file.tellg()) / kRecordBytes;
        const std::int64_t fileDoubles = fileRecords * static_cast<std::int64_t>(kRecordBytes / 8);

        std::int64_t recno = fward;
        std::int64_t hops  = 0;
        while (recno > 0) {
            if (recno >= fileRecords || ++hops > fileRecords) {
                throw Error("SPK: corrupt summary record chain in " + path);
            }
            double summary[kRecordBytes / 8];
            readBytes((recno - 1) * static_cast<std::int64_t>(kRecordBytes), summary, kRecordBytes);
            char names[kRecordBytes];
            readBytes(recno * static_cast<std::int64_t>(kRecordBytes), names, kRecordBytes);
            const double next  = fix(summary[0]);
            const double nsumD = fix(summary[2]);
            if (!(next >= 0.0 && next < static_cast<double>(fileRecords) && next == std::floor(next))
                || !(nsumD >= 0.0 && nsumD <= maxSummaries && nsumD == std::floor(nsumD))) {
                throw Error("SPK: corrupt summary record in " + path);
            }
            const int nsum = static_cast<int>(nsumD);
            for (int i = 0; i < nsum; ++i) {
                const double* s = summary + 3 + i * ss;
                Segment       seg;
                seg.info.startEt = fix(s[0]);
                seg.info.endEt   = fix(s[1]);
                std::int32_t ints[6];
                std::memcpy(ints, s + 2, sizeof(ints));
                for (auto& v : ints) {
                    v = fix(v);
                }
                seg.info.target = ints[0];
                seg.info.center = ints[1];
                seg.info.frame  = ints[2];
                seg.info.type   = ints[3];
                seg.startAddr   = ints[4];
                seg.endAddr     = ints[5];
                std::string name(names + i * ss * 8, static_cast<std::size_t>(ss * 8));
                const auto  last = name.find_last_not_of(std::string(" \0", 2));
                name.erase(last == std::string::npos ? 0 : last + 1);
                seg.info.name = name;
                if (!(std::isfinite(seg.info.startEt) && std::isfinite(seg.info.endEt)
                      && seg.info.startEt <= seg.info.endEt)) {
                    throw Error("SPK: corrupt coverage interval in segment '" + name + "' of " + path);
                }
                if (seg.info.type == 2 || seg.info.type == 3) {
                    readDirectory(seg, fileDoubles);
                }
                infos.push_back(seg.info);
                segs.push_back(std::move(seg));
            }
            recno = static_cast<std::int64_t>(next);
        }
        if (segs.empty()) {
            throw Error("SPK: no segments in " + path);
        }
    }

    // Latest segment for `target` covering `et`.
    const Segment& find(int target, double et) const {
        for (auto it = segs.rbegin(); it != segs.rend(); ++it) {
            if (it->info.target == target && et >= it->info.startEt && et <= it->info.endEt) {
                if (it->info.type != 2 && it->info.type != 3) {
                    throw Error("SPK: unsupported segment type " + std::to_string(it->info.type));
                }
                // States are returned as ICRF axes, which is right for a segment in frame 1 only. One in
                // another frame would need a rotation this reader does not have: refuse it rather than
                // return its numbers as if they were frame 1.
                if (it->info.frame != kFrameJ2000) {
                    throw Error("SPK: segment '" + it->info.name + "' (body " + std::to_string(target)
                                + " relative to body " + std::to_string(it->info.center) + ") of " + path
                                + " is in reference frame " + std::to_string(it->info.frame)
                                + "; only frame 1 (J2000 / ICRF) is supported");
                }
                return *it;
            }
        }
        // Nothing covers `et`: either the epoch is outside the body's coverage, or the kernel does not
        // have the body at all.
        const auto has = [&](int body) {
            return std::any_of(segs.begin(), segs.end(), [&](const Segment& s) { return s.info.target == body; });
        };
        if (has(target)) {
            throw Error("SPK: no segment for body " + std::to_string(target) + " at ET " + std::to_string(et));
        }
        std::string what = "SPK: no segment for body " + std::to_string(target) + " in " + path;
        // A planet (NAIF id n99) and the barycenter of its system (id n) are different bodies, and a
        // kernel may have only the barycenter: DE430 and DE440 have 4 but not 499 (Mars). One is never
        // returned for the other, so say that the barycenter is there to be asked for.
        const int barycenter = target / 100;
        if (target % 100 == 99 && barycenter >= 1 && barycenter <= 9 && has(barycenter)) {
            what += "; the kernel has the barycenter of that planet's system, body " + std::to_string(barycenter)
                  + " (Body::" + std::string(toString(static_cast<Body>(barycenter)))
                  + "), which is a different point: ask for it by name, or load a kernel that has body "
                  + std::to_string(target);
        }
        throw Error(what);
    }

    // Segment state (km, km/s) of seg.target relative to seg.center.
    StateVector evaluate(const Segment& seg, double et) const {
        // Clamp in floating point: casting an out-of-range double to an integer is undefined.
        const double u   = std::floor((et - seg.init) / seg.intlen);
        const auto   idx = static_cast<std::int64_t>(std::clamp(u, 0.0, static_cast<double>(seg.nrec - 1)));
        if (seg.cachedIndex != idx) {
            seg.cache.resize(static_cast<std::size_t>(seg.rsize));
            readDoubles(seg.startAddr + idx * seg.rsize, seg.cache.size(), seg.cache.data());
            seg.cachedIndex = idx;
        }
        const double* r   = seg.cache.data();
        const double  mid = r[0], radius = r[1];
        if (!(std::isfinite(mid) && std::isfinite(radius) && radius > 0.0)) {
            throw Error("SPK: corrupt record in segment '" + seg.info.name + "' of " + path);
        }
        const std::size_t ncomp = seg.info.type == 2 ? 3 : 6;
        const auto        n     = static_cast<std::size_t>((seg.rsize - 2) / static_cast<std::int64_t>(ncomp));
        const double      s     = (et - mid) / radius;
        // Chebyshev polynomials and derivatives.
        std::vector<double> T(n), dT(n);
        T[0]  = 1.0;
        dT[0] = 0.0;
        if (n > 1) {
            T[1]  = s;
            dT[1] = 1.0;
        }
        for (std::size_t k = 2; k < n; ++k) {
            T[k]  = 2.0 * s * T[k - 1] - T[k - 2];
            dT[k] = 2.0 * T[k - 1] + 2.0 * s * dT[k - 1] - dT[k - 2];
        }
        StateVector out;
        for (std::size_t c = 0; c < 3; ++c) {
            const double* coef = r + 2 + c * n;
            double        p = 0.0, v = 0.0;
            for (std::size_t k = n; k-- > 0;) {
                p += coef[k] * T[k];
                v += coef[k] * dT[k];
            }
            out.r[c] = p;
            out.v[c] = v / radius;
            if (ncomp == 6) {
                const double* vc = r + 2 + (3 + c) * n;
                double        vv = 0.0;
                for (std::size_t k = n; k-- > 0;) {
                    vv += vc[k] * T[k];
                }
                out.v[c] = vv;
            }
        }
        return out;
    }

    // Path from body to the root (SSB) as the list of bodies.
    std::vector<int> chain(int body, double et) const {
        std::vector<int> p{body};
        while (body != 0) {
            if (p.size() > kMaxChain) {
                throw Error("SPK: center chain of body " + std::to_string(p.front()) + " does not reach the SSB");
            }
            body = find(body, et).info.center;
            p.push_back(body);
        }
        return p;
    }

    StateVector sumTo(int body, int stop, double et) const {
        StateVector acc;
        for (std::size_t hops = 0; body != stop; ++hops) {
            if (hops > kMaxChain) {
                throw Error("SPK: center chain of body " + std::to_string(body) + " does not reach the SSB");
            }
            const Segment&    seg  = find(body, et);
            const StateVector s    = evaluate(seg, et);
            acc.r                 += s.r;
            acc.v                 += s.v;
            body                   = seg.info.center;
        }
        return acc;
    }
};

SpkEphemeris::SpkEphemeris(const std::string& path)
    : impl_(std::make_unique<Impl>()) {
    impl_->path = path;
    impl_->open();
}

SpkEphemeris::~SpkEphemeris()                                  = default;
SpkEphemeris::SpkEphemeris(SpkEphemeris&&) noexcept            = default;
SpkEphemeris& SpkEphemeris::operator=(SpkEphemeris&&) noexcept = default;

const std::vector<SpkSegmentInfo>& SpkEphemeris::segments() const noexcept {
    return impl_->infos;
}
bool SpkEphemeris::bigEndian() const noexcept {
    return impl_->big;
}

StateVector SpkEphemeris::stateById(int target, int center, double et) const {
    std::lock_guard lock(impl_->mutex);
    if (target == center) {
        return {};
    }
    // Find the first common node of the two paths to the SSB.
    const std::vector<int> pt     = target == 0 ? std::vector<int>{0} : impl_->chain(target, et);
    const std::vector<int> pc     = center == 0 ? std::vector<int>{0} : impl_->chain(center, et);
    int                    common = 0;
    for (int b : pt) {
        if (std::find(pc.begin(), pc.end(), b) != pc.end()) {
            common = b;
            break;
        }
    }
    const StateVector a = impl_->sumTo(target, common, et);
    const StateVector b = impl_->sumTo(center, common, et);
    return {(a.r - b.r) * 1000.0, (a.v - b.v) * 1000.0};
}

StateVector SpkEphemeris::state(Body target, Body center, const Time& t) const {
    return stateById(naifId(target), naifId(center), ephemerisTime(t));
}

double ephemerisTime(const Time& t) {
    return t.to(TimeScale::TDB).jd().secondsSince(JulianDate{constants::jdJ2000, 0.0});
}

} // namespace astro
