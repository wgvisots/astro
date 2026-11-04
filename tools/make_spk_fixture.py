#!/usr/bin/env python3
###########################################################################
# make_spk_fixture.py
#
# Copyright 2024 W. Greg Visotski
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
###########################################################################
"""Cut a small SPK test fixture out of a JPL planetary ephemeris.

Usage: python3 tools/make_spk_fixture.py path/to/de440s.bsp

Writes tests/data/de440s_2024q1_le.bsp (little-endian) and _be.bsp (big-endian)
containing 2024-01-01..2024-03-01 for segments 0->10, 0->3, 3->301, 3->399 as
type 2, plus 0->4 re-encoded as type 3 (position + velocity Chebyshev), and
tests/data/spk_reference.csv with states computed by jplephem from the full
kernel. Requires numpy and jplephem (dev only). JPL ephemerides are public domain.
"""
import os
import struct
import sys

import numpy as np
from numpy.polynomial import chebyshev as C
from jplephem.spk import SPK

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tests", "data")
J2000 = 2451545.0
START_JD, END_JD = 2460310.5, 2460371.5  # 2024-01-01 .. 2024-03-01 TDB
FTPSTR = b"FTPSTR:\r:\n:\r\n:\r\x00:\x81:\x10\xce:ENDFTP"


def subset_records(seg):
    """Return (init, intlen, deg, records[n, 3, deg+1]) covering the window."""
    init_jd, intlen_days, coef = seg.load_array()  # coef: (3, N, deg+1); jplephem returns JD and days
    init = (init_jd - J2000) * 86400.0
    intlen = intlen_days * 86400.0
    t0 = (START_JD - J2000) * 86400.0
    t1 = (END_JD - J2000) * 86400.0
    i0 = int((t0 - init) // intlen)
    i1 = int(np.ceil((t1 - init) / intlen))
    recs = np.transpose(coef[:, i0:i1, :], (1, 0, 2))  # (n, 3, deg+1)
    return init + i0 * intlen, intlen, coef.shape[2] - 1, recs


def segment_data(seg, as_type3=False):
    init, intlen, deg, recs = subset_records(seg)
    rows = []
    for k, rec in enumerate(recs):
        mid = init + (k + 0.5) * intlen
        radius = intlen / 2.0
        row = [mid, radius] + list(rec[0]) + list(rec[1]) + list(rec[2])
        if as_type3:
            for comp in rec:
                d = C.chebder(comp) / radius  # km/s per Chebyshev coefficient
                row += list(np.concatenate([d, [0.0]]))
        rows.append(row)
    rsize = len(rows[0])
    data = [x for row in rows for x in row] + [init, intlen, float(rsize), float(len(rows))]
    start_et = init
    end_et = init + len(rows) * intlen
    return data, start_et, end_et


def write_daf(path, segments, big_endian):
    e = ">" if big_endian else "<"
    nd, ni = 2, 6
    ss = nd + (ni + 1) // 2
    # Records: 1 file record, 2 summary, 3 names, data from record 4 (address 3*128+1).
    addr = 3 * 128 + 1
    summaries, data_words = [], []
    for (target, center, typ, data, start_et, end_et) in segments:
        start = addr
        end = addr + len(data) - 1
        summaries.append((start_et, end_et, target, center, 1, typ, start, end))
        data_words += data
        addr = end + 1
    free = addr
    fr = bytearray(1024)
    fr[0:8] = b"DAF/SPK "
    fr[8:16] = struct.pack(e + "ii", nd, ni)
    fr[16:76] = b"astro test fixture (DE440s subset)".ljust(60)
    fr[76:88] = struct.pack(e + "iii", 2, 2, free)
    fr[88:96] = b"BIG-IEEE" if big_endian else b"LTL-IEEE"
    fr[699:699 + len(FTPSTR)] = FTPSTR
    sr = bytearray(1024)
    sr[0:24] = struct.pack(e + "ddd", 0.0, 0.0, float(len(summaries)))
    off = 24
    for s in summaries:
        sr[off:off + 8 * ss] = struct.pack(e + "dd6i", *s)
        off += 8 * ss
    nr = bytearray(b" " * 1024)
    for i, s in enumerate(summaries):
        name = f"SEG {s[3]}->{s[2]} type {s[5]}".encode().ljust(8 * ss)
        nr[i * 8 * ss:(i + 1) * 8 * ss] = name
    body = struct.pack(e + f"{len(data_words)}d", *data_words)
    body += b"\0" * ((-len(body)) % 1024)
    with open(path, "wb") as f:
        f.write(bytes(fr) + bytes(sr) + bytes(nr) + body)
    print(f"wrote {path} ({os.path.getsize(path)} bytes)")


def main():
    full = SPK.open(sys.argv[1])
    wanted = [(0, 10, False), (0, 3, False), (3, 301, False), (3, 399, False), (0, 4, True)]
    segs = []
    for center, target, t3 in wanted:
        seg = full[center, target]
        data, s0, s1 = segment_data(seg, as_type3=t3)
        segs.append((target, center, 3 if t3 else 2, data, s0, s1))
    le = os.path.join(OUT, "de440s_2024q1_le.bsp")
    be = os.path.join(OUT, "de440s_2024q1_be.bsp")
    write_daf(le, segs, False)
    write_daf(be, segs, True)

    # Verify the fixture with jplephem and write reference states from the full kernel.
    fix = SPK.open(le)
    rows = []
    for jd in np.linspace(START_JD + 0.1, END_JD - 0.1, 25):
        for center, target, _ in wanted:
            p, v = full[center, target].compute_and_differentiate(jd)
            fseg = fix[center, target]
            if fseg.data_type == 3:  # type 3 stores velocity (km/s) explicitly
                pv = fseg.compute(jd)
                pf, vf = pv[:3], pv[3:] * 86400.0
            else:
                pf, vf = fseg.compute_and_differentiate(jd)
            # ~1e-13 relative differences arise from epoch rounding in the re-encoded records.
            assert np.max(np.abs(p - pf)) < 1e-12 * np.linalg.norm(p) + 1e-6, (center, target, p - pf)
            assert np.max(np.abs(v - vf)) < 1e-6, (center, target, v - vf)
            rows.append([jd, center, target, *p, *(v / 86400.0)])
    with open(os.path.join(OUT, "spk_reference.csv"), "w") as f:
        f.write("# DE440s states from jplephem (km, km/s, ICRF), TDB Julian dates\n")
        f.write("jd_tdb,center,target,x,y,z,vx,vy,vz\n")
        for r in rows:
            f.write(",".join([repr(float(r[0])), str(r[1]), str(r[2])] + [repr(float(x)) for x in r[3:]]) + "\n")
    print(f"wrote spk_reference.csv ({len(rows)} rows); fixture verified with jplephem")


if __name__ == "__main__":
    main()
