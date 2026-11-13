/**********************************************************************
 * sgp4.cpp
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

// SGP4/SDP4, ported from the reference C++ implementation published with Vallado, Crawford, Hujsak &
// Kelso (2006), "Revisiting Spacetrack Report #3", AIAA 2006-6753
// (https://celestrak.org/publications/AIAA/2006-6753/). The routines (dpper, dscom, dsinit, dspace,
// initl, sgp4) and their variable names are kept so the code can be compared with the reference line
// by line, and the published "sgp4fix" corrections are included. See THIRD_PARTY_NOTICES.md.

#include "astro/orbit/sgp4.h"

#include <cmath>
#include <string>

#include "../core/enum_io.h"

namespace astro {

namespace {
constexpr double pi      = 3.14159265358979323846;
constexpr double twopi   = 2.0 * pi;
constexpr double x2o3    = 2.0 / 3.0;
constexpr double deg2rad = pi / 180.0;
constexpr double xpdotp  = 1440.0 / twopi; // rev/day -> rad/min

double gstime(double jdut1) {
    const double tut1 = (jdut1 - 2451545.0) / 36525.0;
    double temp = -6.2e-6 * tut1 * tut1 * tut1 + 0.093104 * tut1 * tut1 + (876600.0 * 3600.0 + 8640184.812866) * tut1
                + 67310.54841; // seconds
    temp = std::fmod(temp * deg2rad / 240.0, twopi);
    if (temp < 0.0) {
        temp += twopi;
    }
    return temp;
}

const char* errorMessage(int code) {
    switch (code) {
    case 1:  return "SGP4: mean eccentricity out of range (e >= 1 or e < -0.001)";
    case 2:  return "SGP4: mean motion less than zero";
    case 3:  return "SGP4: perturbed eccentricity out of range";
    case 4:  return "SGP4: semi-latus rectum less than zero";
    case 6:  return "SGP4: satellite has decayed";
    default: return "SGP4: propagation error";
    }
}

// Throws astro::Error naming the element if it is NaN or infinite. The range tests of the reference code
// (e >= 1, n <= 0, ...) are all false for NaN, so a non-finite element would pass them and come out as a NaN state
// with error code 0.
void requireFinite(double value, const char* name) {
    if (!std::isfinite(value)) {
        throw Error(std::string("Sgp4: ") + name + " is not finite");
    }
}

// (86400 s) / (length of the UTC day that starts at the Julian date `midnight`): the factor between the day fraction
// of astro::Time, which is the time of day over the length of that day (86401 s when it ends with a leap second),
// and the day fraction of an element set, which is always of 86400 s. 12:00:00 is 43200 s into any day, so the factor
// is twice the fraction at which astro::Time puts noon. It is exactly 1 on an ordinary day. (A copy of the function
// in tle.cpp, which reads and writes the epoch field with it.)
double tleDayPerUtcDay(double midnight) {
    int year = 0, month = 0, day = 0;
    dateFromJulianDayNumber(std::llround(midnight + 0.5), year, month, day);
    return 2.0 * Time::fromCalendar(TimeScale::UTC, year, month, day, 12).jd().jd2;
}

// Outputs of dscom needed later.
struct DsCom {
    double snodm, cnodm, sinim, cosim, sinomm, cosomm, day, e3, ee2, em, emsq, gam, peo, pgho, pho, pinco, plo, rtemsq,
        se2, se3, sgh2, sgh3, sgh4, sh2, sh3, si2, si3, sl2, sl3, sl4, s1, s2, s3, s4, s5, s6, s7, ss1, ss2, ss3, ss4,
        ss5, ss6, ss7, sz1, sz2, sz3, sz11, sz12, sz13, sz21, sz22, sz23, sz31, sz32, sz33, xgh2, xgh3, xgh4, xh2, xh3,
        xi2, xi3, xl2, xl3, xl4, nm, z1, z2, z3, z11, z12, z13, z21, z22, z23, z31, z32, z33, zmol, zmos;
};

} // namespace

struct Sgp4::Impl {
    Tle  tle;
    Time epochTime;
    char opsmode{'i'};

    // Gravity constants.
    double tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2;

    // Elements.
    double bstar, ndot, nddot, ecco, argpo, inclo, mo, no_kozai, nodeo, no_unkozai, epochDays;

    // Near-earth.
    int    isimp{0};
    char   method{'n'};
    double aycof{0}, con41{0}, cc1{0}, cc4{0}, cc5{0}, d2{0}, d3{0}, d4{0}, delmo{0}, eta{0}, argpdot{0}, omgcof{0},
        sinmao{0}, t2cof{0}, t3cof{0}, t4cof{0}, t5cof{0}, x1mth2{0}, x7thm1{0}, mdot{0}, nodedot{0}, xlcof{0},
        xmcof{0}, nodecf{0};

    // Deep space.
    int    irez{0};
    double d2201{0}, d2211{0}, d3210{0}, d3222{0}, d4410{0}, d4422{0}, d5220{0}, d5232{0}, d5421{0}, d5433{0}, dedt{0},
        del1{0}, del2{0}, del3{0}, didt{0}, dmdt{0}, dnodt{0}, domdt{0}, e3{0}, ee2{0}, peo{0}, pgho{0}, pho{0},
        pinco{0}, plo{0}, se2{0}, se3{0}, sgh2{0}, sgh3{0}, sgh4{0}, sh2{0}, sh3{0}, si2{0}, si3{0}, sl2{0}, sl3{0},
        sl4{0}, gsto{0}, xfact{0}, xgh2{0}, xgh3{0}, xgh4{0}, xh2{0}, xh3{0}, xi2{0}, xi3{0}, xl2{0}, xl3{0}, xl4{0},
        xlamo{0}, zmol{0}, zmos{0}, atime{0}, xli{0}, xni{0};

    void setGravity(Sgp4Gravity g) {
        switch (g) {
        case Sgp4Gravity::WGS72Old:
            mu            = 398600.79964;
            radiusearthkm = 6378.135;
            xke           = 0.0743669161;
            tumin         = 1.0 / xke;
            j2            = 0.001082616;
            j3            = -0.00000253881;
            j4            = -0.00000165597;
            break;
        case Sgp4Gravity::WGS72:
            mu            = 398600.8;
            radiusearthkm = 6378.135;
            xke           = 60.0 / std::sqrt(radiusearthkm * radiusearthkm * radiusearthkm / mu);
            tumin         = 1.0 / xke;
            j2            = 0.001082616;
            j3            = -0.00000253881;
            j4            = -0.00000165597;
            break;
        case Sgp4Gravity::WGS84:
            mu            = 398600.5;
            radiusearthkm = 6378.137;
            xke           = 60.0 / std::sqrt(radiusearthkm * radiusearthkm * radiusearthkm / mu);
            tumin         = 1.0 / xke;
            j2            = 0.00108262998905;
            j3            = -0.00000253215306;
            j4            = -0.00000161098761;
            break;
        }
        j3oj2 = j3 / j2;
    }

    // ----------------------------------------------------------- dpper
    void dpper(double t, bool init, double& ep, double& inclp, double& nodep, double& argpp, double& mp) const {
        constexpr double zns = 1.19459e-5, zes = 0.01675, znl = 1.5835218e-4, zel = 0.05490;
        double           zm    = init ? zmos : zmos + zns * t;
        double           zf    = zm + 2.0 * zes * std::sin(zm);
        double           sinzf = std::sin(zf);
        double           f2    = 0.5 * sinzf * sinzf - 0.25;
        double           f3    = -0.5 * sinzf * std::cos(zf);
        const double     ses   = se2 * f2 + se3 * f3;
        const double     sis   = si2 * f2 + si3 * f3;
        const double     sls   = sl2 * f2 + sl3 * f3 + sl4 * sinzf;
        const double     sghs  = sgh2 * f2 + sgh3 * f3 + sgh4 * sinzf;
        const double     shs   = sh2 * f2 + sh3 * f3;
        zm                     = init ? zmol : zmol + znl * t;
        zf                     = zm + 2.0 * zel * std::sin(zm);
        sinzf                  = std::sin(zf);
        f2                     = 0.5 * sinzf * sinzf - 0.25;
        f3                     = -0.5 * sinzf * std::cos(zf);
        const double sel       = ee2 * f2 + e3 * f3;
        const double sil       = xi2 * f2 + xi3 * f3;
        const double sll       = xl2 * f2 + xl3 * f3 + xl4 * sinzf;
        const double sghl      = xgh2 * f2 + xgh3 * f3 + xgh4 * sinzf;
        const double shll      = xh2 * f2 + xh3 * f3;
        double       pe        = ses + sel;
        double       pinc      = sis + sil;
        double       pl        = sls + sll;
        double       pgh       = sghs + sghl;
        double       ph        = shs + shll;
        if (init) {
            return;
        }

        pe                 -= peo;
        pinc               -= pinco;
        pl                 -= plo;
        pgh                -= pgho;
        ph                 -= pho;
        inclp              += pinc;
        ep                 += pe;
        const double sinip = std::sin(inclp), cosip = std::cos(inclp);
        // sgp4fix: use the perturbed inclination for the Lyddane choice (GSFC approach).
        if (inclp >= 0.2) {
            ph    /= sinip;
            pgh   -= cosip * ph;
            argpp += pgh;
            nodep += ph;
            mp    += pl;
        } else {
            // Lyddane modification.
            const double sinop = std::sin(nodep), cosop = std::cos(nodep);
            double       alfdp  = sinip * sinop;
            double       betdp  = sinip * cosop;
            const double dalf   = ph * cosop + pinc * cosip * sinop;
            const double dbet   = -ph * sinop + pinc * cosip * cosop;
            alfdp              += dalf;
            betdp              += dbet;
            nodep               = std::fmod(nodep, twopi);
            if (nodep < 0.0 && opsmode == 'a') {
                nodep += twopi;
            }
            double       xls   = mp + argpp + cosip * nodep;
            const double dls   = pl + pgh - pinc * nodep * sinip;
            xls               += dls;
            const double xnoh  = nodep;
            nodep              = std::atan2(alfdp, betdp);
            if (nodep < 0.0 && opsmode == 'a') {
                nodep += twopi;
            }
            if (std::abs(xnoh - nodep) > pi) {
                if (nodep < xnoh) {
                    nodep += twopi;
                } else {
                    nodep -= twopi;
                }
            }
            mp    += pl;
            argpp  = xls - mp - cosip * nodep;
        }
    }

    // ----------------------------------------------------------- dscom
    static DsCom dscom(double epoch, double ep, double argpp, double tc, double inclp, double nodep, double np) {
        constexpr double zes = 0.01675, zel = 0.05490, c1ss = 2.9864797e-6, c1l = 4.7968065e-7, zsinis = 0.39785416,
                         zcosis = 0.91744867, zcosgs = 0.1945905, zsings = -0.98088458;
        DsCom d{};
        d.nm                = np;
        d.em                = ep;
        d.snodm             = std::sin(nodep);
        d.cnodm             = std::cos(nodep);
        d.sinomm            = std::sin(argpp);
        d.cosomm            = std::cos(argpp);
        d.sinim             = std::sin(inclp);
        d.cosim             = std::cos(inclp);
        d.emsq              = d.em * d.em;
        const double betasq = 1.0 - d.emsq;
        d.rtemsq            = std::sqrt(betasq);

        d.peo = d.pinco = d.plo = d.pgho = d.pho = 0.0;
        d.day                                    = epoch + 18261.5 + tc / 1440.0;
        const double xnodce                      = std::fmod(4.5236020 - 9.2422029e-4 * d.day, twopi);
        const double stem = std::sin(xnodce), ctem = std::cos(xnodce);
        const double zcosil = 0.91375164 - 0.03568096 * ctem;
        const double zsinil = std::sqrt(1.0 - zcosil * zcosil);
        const double zsinhl = 0.089683511 * stem / zsinil;
        const double zcoshl = std::sqrt(1.0 - zsinhl * zsinhl);
        d.gam               = 5.8351514 + 0.0019443680 * d.day;
        double       zx     = 0.39785416 * stem / zsinil;
        const double zy     = zcoshl * ctem + 0.91744867 * zsinhl * stem;
        zx                  = std::atan2(zx, zy);
        zx                  = d.gam + zx - xnodce;
        const double zcosgl = std::cos(zx), zsingl = std::sin(zx);

        double       zcosg = zcosgs, zsing = zsings, zcosi = zcosis, zsini = zsinis, zcosh = d.cnodm, zsinh = d.snodm;
        double       cc   = c1ss;
        const double xnoi = 1.0 / d.nm;
        for (int lsflg = 1; lsflg <= 2; ++lsflg) {
            const double a1  = zcosg * zcosh + zsing * zcosi * zsinh;
            const double a3  = -zsing * zcosh + zcosg * zcosi * zsinh;
            const double a7  = -zcosg * zsinh + zsing * zcosi * zcosh;
            const double a8  = zsing * zsini;
            const double a9  = zsing * zsinh + zcosg * zcosi * zcosh;
            const double a10 = zcosg * zsini;
            const double a2  = d.cosim * a7 + d.sinim * a8;
            const double a4  = d.cosim * a9 + d.sinim * a10;
            const double a5  = -d.sinim * a7 + d.cosim * a8;
            const double a6  = -d.sinim * a9 + d.cosim * a10;
            const double x1  = a1 * d.cosomm + a2 * d.sinomm;
            const double x2  = a3 * d.cosomm + a4 * d.sinomm;
            const double x3  = -a1 * d.sinomm + a2 * d.cosomm;
            const double x4  = -a3 * d.sinomm + a4 * d.cosomm;
            const double x5  = a5 * d.sinomm;
            const double x6  = a6 * d.sinomm;
            const double x7  = a5 * d.cosomm;
            const double x8  = a6 * d.cosomm;
            d.z31            = 12.0 * x1 * x1 - 3.0 * x3 * x3;
            d.z32            = 24.0 * x1 * x2 - 6.0 * x3 * x4;
            d.z33            = 12.0 * x2 * x2 - 3.0 * x4 * x4;
            d.z1             = 3.0 * (a1 * a1 + a2 * a2) + d.z31 * d.emsq;
            d.z2             = 6.0 * (a1 * a3 + a2 * a4) + d.z32 * d.emsq;
            d.z3             = 3.0 * (a3 * a3 + a4 * a4) + d.z33 * d.emsq;
            d.z11            = -6.0 * a1 * a5 + d.emsq * (-24.0 * x1 * x7 - 6.0 * x3 * x5);
            d.z12 = -6.0 * (a1 * a6 + a3 * a5) + d.emsq * (-24.0 * (x2 * x7 + x1 * x8) - 6.0 * (x3 * x6 + x4 * x5));
            d.z13 = -6.0 * a3 * a6 + d.emsq * (-24.0 * x2 * x8 - 6.0 * x4 * x6);
            d.z21 = 6.0 * a2 * a5 + d.emsq * (24.0 * x1 * x5 - 6.0 * x3 * x7);
            d.z22 = 6.0 * (a4 * a5 + a2 * a6) + d.emsq * (24.0 * (x2 * x5 + x1 * x6) - 6.0 * (x4 * x7 + x3 * x8));
            d.z23 = 6.0 * a4 * a6 + d.emsq * (24.0 * x2 * x6 - 6.0 * x4 * x8);
            d.z1  = d.z1 + d.z1 + betasq * d.z31;
            d.z2  = d.z2 + d.z2 + betasq * d.z32;
            d.z3  = d.z3 + d.z3 + betasq * d.z33;
            d.s3  = cc * xnoi;
            d.s2  = -0.5 * d.s3 / d.rtemsq;
            d.s4  = d.s3 * d.rtemsq;
            d.s1  = -15.0 * d.em * d.s4;
            d.s5  = x1 * x3 + x2 * x4;
            d.s6  = x2 * x3 + x1 * x4;
            d.s7  = x2 * x4 - x1 * x3;
            if (lsflg == 1) {
                d.ss1  = d.s1;
                d.ss2  = d.s2;
                d.ss3  = d.s3;
                d.ss4  = d.s4;
                d.ss5  = d.s5;
                d.ss6  = d.s6;
                d.ss7  = d.s7;
                d.sz1  = d.z1;
                d.sz2  = d.z2;
                d.sz3  = d.z3;
                d.sz11 = d.z11;
                d.sz12 = d.z12;
                d.sz13 = d.z13;
                d.sz21 = d.z21;
                d.sz22 = d.z22;
                d.sz23 = d.z23;
                d.sz31 = d.z31;
                d.sz32 = d.z32;
                d.sz33 = d.z33;
                zcosg  = zcosgl;
                zsing  = zsingl;
                zcosi  = zcosil;
                zsini  = zsinil;
                zcosh  = zcoshl * d.cnodm + zsinhl * d.snodm;
                zsinh  = d.snodm * zcoshl - d.cnodm * zsinhl;
                cc     = c1l;
            }
        }
        d.zmol = std::fmod(4.7199672 + 0.22997150 * d.day - d.gam, twopi);
        d.zmos = std::fmod(6.2565837 + 0.017201977 * d.day, twopi);

        // Solar terms.
        d.se2  = 2.0 * d.ss1 * d.ss6;
        d.se3  = 2.0 * d.ss1 * d.ss7;
        d.si2  = 2.0 * d.ss2 * d.sz12;
        d.si3  = 2.0 * d.ss2 * (d.sz13 - d.sz11);
        d.sl2  = -2.0 * d.ss3 * d.sz2;
        d.sl3  = -2.0 * d.ss3 * (d.sz3 - d.sz1);
        d.sl4  = -2.0 * d.ss3 * (-21.0 - 9.0 * d.emsq) * zes;
        d.sgh2 = 2.0 * d.ss4 * d.sz32;
        d.sgh3 = 2.0 * d.ss4 * (d.sz33 - d.sz31);
        d.sgh4 = -18.0 * d.ss4 * zes;
        d.sh2  = -2.0 * d.ss2 * d.sz22;
        d.sh3  = -2.0 * d.ss2 * (d.sz23 - d.sz21);
        // Lunar terms.
        d.ee2  = 2.0 * d.s1 * d.s6;
        d.e3   = 2.0 * d.s1 * d.s7;
        d.xi2  = 2.0 * d.s2 * d.z12;
        d.xi3  = 2.0 * d.s2 * (d.z13 - d.z11);
        d.xl2  = -2.0 * d.s3 * d.z2;
        d.xl3  = -2.0 * d.s3 * (d.z3 - d.z1);
        d.xl4  = -2.0 * d.s3 * (-21.0 - 9.0 * d.emsq) * zel;
        d.xgh2 = 2.0 * d.s4 * d.z32;
        d.xgh3 = 2.0 * d.s4 * (d.z33 - d.z31);
        d.xgh4 = -18.0 * d.s4 * zel;
        d.xh2  = -2.0 * d.s2 * d.z22;
        d.xh3  = -2.0 * d.s2 * (d.z23 - d.z21);
        return d;
    }

    // ----------------------------------------------------------- dsinit
    void dsinit(const DsCom& d,
                double       t,
                double       tc,
                double       xpidot,
                double&      em,
                double&      argpm,
                double&      inclm,
                double&      mm,
                double&      nm,
                double&      nodem,
                double&      dndt) {
        constexpr double q22 = 1.7891679e-6, q31 = 2.1460748e-6, q33 = 2.2123015e-7, root22 = 1.7891679e-6,
                         root44 = 7.3636953e-9, root54 = 2.1765803e-9, rptim = 4.37526908801129966e-3,
                         root32 = 3.7393792e-7, root52 = 1.1428639e-7, znl = 1.5835218e-4, zns = 1.19459e-5;
        const double cosim = d.cosim, sinim = d.sinim;
        double       emsq = d.emsq;

        irez = 0;
        if (nm < 0.0052359877 && nm > 0.0034906585) {
            irez = 1;
        }
        if (nm >= 8.26e-3 && nm <= 9.24e-3 && em >= 0.5) {
            irez = 2;
        }

        // Solar terms.
        const double ses  = d.ss1 * zns * d.ss5;
        const double sis  = d.ss2 * zns * (d.sz11 + d.sz13);
        const double sls  = -zns * d.ss3 * (d.sz1 + d.sz3 - 14.0 - 6.0 * emsq);
        const double sghs = d.ss4 * zns * (d.sz31 + d.sz33 - 6.0);
        double       shs  = -zns * d.ss2 * (d.sz21 + d.sz23);
        if (inclm < 5.2359877e-2 || inclm > pi - 5.2359877e-2) {
            shs = 0.0;
        }
        if (sinim != 0.0) {
            shs /= sinim;
        }
        const double sgs = sghs - cosim * shs;

        // Lunar terms.
        dedt              = ses + d.s1 * znl * d.s5;
        didt              = sis + d.s2 * znl * (d.z11 + d.z13);
        dmdt              = sls - znl * d.s3 * (d.z1 + d.z3 - 14.0 - 6.0 * emsq);
        const double sghl = d.s4 * znl * (d.z31 + d.z33 - 6.0);
        double       shll = -znl * d.s2 * (d.z21 + d.z23);
        if (inclm < 5.2359877e-2 || inclm > pi - 5.2359877e-2) {
            shll = 0.0;
        }
        domdt = sgs + sghl;
        dnodt = shs;
        if (sinim != 0.0) {
            domdt -= cosim / sinim * shll;
            dnodt += shll / sinim;
        }

        // Deep-space resonance effects.
        dndt                = 0.0;
        const double theta  = std::fmod(gsto + tc * rptim, twopi);
        em                 += dedt * t;
        inclm              += didt * t;
        argpm              += domdt * t;
        nodem              += dnodt * t;
        mm                 += dmdt * t;

        if (irez != 0) {
            const double aonv = std::pow(nm / xke, x2o3);
            if (irez == 2) {
                // Geopotential resonance for 12-hour orbits.
                const double cosisq = cosim * cosim;
                const double emo    = em;
                em                  = ecco;
                const double emsqo  = emsq;
                emsq                = ecco * ecco;
                const double eoc    = em * emsq;
                const double g201   = -0.306 - (em - 0.64) * 0.440;
                double       g211, g310, g322, g410, g422, g520, g521, g532, g533;
                if (em <= 0.65) {
                    g211 = 3.616 - 13.2470 * em + 16.2900 * emsq;
                    g310 = -19.302 + 117.3900 * em - 228.4190 * emsq + 156.5910 * eoc;
                    g322 = -18.9068 + 109.7927 * em - 214.6334 * emsq + 146.5816 * eoc;
                    g410 = -41.122 + 242.6940 * em - 471.0940 * emsq + 313.9530 * eoc;
                    g422 = -146.407 + 841.8800 * em - 1629.014 * emsq + 1083.4350 * eoc;
                    g520 = -532.114 + 3017.977 * em - 5740.032 * emsq + 3708.2760 * eoc;
                } else {
                    g211 = -72.099 + 331.819 * em - 508.738 * emsq + 266.724 * eoc;
                    g310 = -346.844 + 1582.851 * em - 2415.925 * emsq + 1246.113 * eoc;
                    g322 = -342.585 + 1554.908 * em - 2366.899 * emsq + 1215.972 * eoc;
                    g410 = -1052.797 + 4758.686 * em - 7193.992 * emsq + 3651.957 * eoc;
                    g422 = -3581.690 + 16178.110 * em - 24462.770 * emsq + 12422.520 * eoc;
                    if (em > 0.715) {
                        g520 = -5149.66 + 29936.92 * em - 54087.36 * emsq + 31324.56 * eoc;
                    } else {
                        g520 = 1464.74 - 4664.75 * em + 3763.64 * emsq;
                    }
                }
                if (em < 0.7) {
                    g533 = -919.22770 + 4988.6100 * em - 9064.7700 * emsq + 5542.21 * eoc;
                    g521 = -822.71072 + 4568.6173 * em - 8491.4146 * emsq + 5337.524 * eoc;
                    g532 = -853.66600 + 4690.2500 * em - 8624.7700 * emsq + 5341.4 * eoc;
                } else {
                    g533 = -37995.780 + 161616.52 * em - 229838.20 * emsq + 109377.94 * eoc;
                    g521 = -51752.104 + 218913.95 * em - 309468.16 * emsq + 146349.42 * eoc;
                    g532 = -40023.880 + 170470.89 * em - 242699.48 * emsq + 115605.82 * eoc;
                }
                const double sini2 = sinim * sinim;
                const double f220  = 0.75 * (1.0 + 2.0 * cosim + cosisq);
                const double f221  = 1.5 * sini2;
                const double f321  = 1.875 * sinim * (1.0 - 2.0 * cosim - 3.0 * cosisq);
                const double f322  = -1.875 * sinim * (1.0 + 2.0 * cosim - 3.0 * cosisq);
                const double f441  = 35.0 * sini2 * f220;
                const double f442  = 39.3750 * sini2 * sini2;
                const double f522 =
                    9.84375 * sinim
                    * (sini2 * (1.0 - 2.0 * cosim - 5.0 * cosisq) + 0.33333333 * (-2.0 + 4.0 * cosim + 6.0 * cosisq));
                const double f523 = sinim
                                  * (4.92187512 * sini2 * (-2.0 - 4.0 * cosim + 10.0 * cosisq)
                                     + 6.56250012 * (1.0 + 2.0 * cosim - 3.0 * cosisq));
                const double f542 =
                    29.53125 * sinim * (2.0 - 8.0 * cosim + cosisq * (-12.0 + 8.0 * cosim + 10.0 * cosisq));
                const double f543 =
                    29.53125 * sinim * (-2.0 - 8.0 * cosim + cosisq * (12.0 + 8.0 * cosim - 10.0 * cosisq));
                const double xno2   = nm * nm;
                const double ainv2  = aonv * aonv;
                double       temp1  = 3.0 * xno2 * ainv2;
                double       temp   = temp1 * root22;
                d2201               = temp * f220 * g201;
                d2211               = temp * f221 * g211;
                temp1              *= aonv;
                temp                = temp1 * root32;
                d3210               = temp * f321 * g310;
                d3222               = temp * f322 * g322;
                temp1              *= aonv;
                temp                = 2.0 * temp1 * root44;
                d4410               = temp * f441 * g410;
                d4422               = temp * f442 * g422;
                temp1              *= aonv;
                temp                = temp1 * root52;
                d5220               = temp * f522 * g520;
                d5232               = temp * f523 * g532;
                temp                = 2.0 * temp1 * root54;
                d5421               = temp * f542 * g521;
                d5433               = temp * f543 * g533;
                xlamo               = std::fmod(mo + nodeo + nodeo - theta - theta, twopi);
                xfact               = mdot + dmdt + 2.0 * (nodedot + dnodt - rptim) - no_unkozai;
                em                  = emo;
                emsq                = emsqo;
            }
            if (irez == 1) {
                // Synchronous resonance terms.
                const double g200 = 1.0 + emsq * (-2.5 + 0.8125 * emsq);
                const double g310 = 1.0 + 2.0 * emsq;
                const double g300 = 1.0 + emsq * (-6.0 + 6.60937 * emsq);
                const double f220 = 0.75 * (1.0 + cosim) * (1.0 + cosim);
                const double f311 = 0.9375 * sinim * sinim * (1.0 + 3.0 * cosim) - 0.75 * (1.0 + cosim);
                double       f330 = 1.0 + cosim;
                f330              = 1.875 * f330 * f330 * f330;
                del1              = 3.0 * nm * nm * aonv * aonv;
                del2              = 2.0 * del1 * f220 * g200 * q22;
                del3              = 3.0 * del1 * f330 * g300 * q33 * aonv;
                del1              = del1 * f311 * g310 * q31 * aonv;
                xlamo             = std::fmod(mo + nodeo + argpo - theta, twopi);
                xfact             = mdot + xpidot - rptim + dmdt + domdt + dnodt - no_unkozai;
            }
            // Initialize the integrator.
            xli   = xlamo;
            xni   = no_unkozai;
            atime = 0.0;
            nm    = no_unkozai + dndt;
        }
    }

    // ----------------------------------------------------------- dspace
    void dspace(double  t,
                double  tc,
                double& atimeL,
                double& em,
                double& argpm,
                double& inclm,
                double& xliL,
                double& mm,
                double& xniL,
                double& nodem,
                double& dndt,
                double& nm) const {
        constexpr double fasx2 = 0.13130908, fasx4 = 2.8843198, fasx6 = 0.37448087, g22 = 5.7686396, g32 = 0.95240898,
                         g44 = 1.8014998, g52 = 1.0508330, g54 = 4.4108898, rptim = 4.37526908801129966e-3,
                         stepp = 720.0, stepn = -720.0, step2 = 259200.0;
        dndt                = 0.0;
        const double theta  = std::fmod(gsto + tc * rptim, twopi);
        em                 += dedt * t;
        inclm              += didt * t;
        argpm              += domdt * t;
        nodem              += dnodt * t;
        mm                 += dmdt * t;

        double ft = 0.0;
        if (irez != 0) {
            if (atimeL == 0.0 || t * atimeL <= 0.0 || std::abs(t) < std::abs(atimeL)) {
                atimeL = 0.0;
                xniL   = no_unkozai;
                xliL   = xlamo;
            }
            const double delt = t > 0.0 ? stepp : stepn;
            double       xndt = 0.0, xldot = 0.0, xnddt = 0.0;
            for (;;) {
                if (irez != 2) {
                    // Near-synchronous resonance terms.
                    xndt = del1 * std::sin(xliL - fasx2) + del2 * std::sin(2.0 * (xliL - fasx4))
                         + del3 * std::sin(3.0 * (xliL - fasx6));
                    xldot = xniL + xfact;
                    xnddt = del1 * std::cos(xliL - fasx2) + 2.0 * del2 * std::cos(2.0 * (xliL - fasx4))
                          + 3.0 * del3 * std::cos(3.0 * (xliL - fasx6));
                    xnddt *= xldot;
                } else {
                    // Near-half-day resonance terms.
                    const double xomi  = argpo + argpdot * atimeL;
                    const double x2omi = xomi + xomi;
                    const double x2li  = xliL + xliL;
                    xndt               = d2201 * std::sin(x2omi + xliL - g22) + d2211 * std::sin(xliL - g22)
                         + d3210 * std::sin(xomi + xliL - g32) + d3222 * std::sin(-xomi + xliL - g32)
                         + d4410 * std::sin(x2omi + x2li - g44) + d4422 * std::sin(x2li - g44)
                         + d5220 * std::sin(xomi + xliL - g52) + d5232 * std::sin(-xomi + xliL - g52)
                         + d5421 * std::sin(xomi + x2li - g54) + d5433 * std::sin(-xomi + x2li - g54);
                    xldot = xniL + xfact;
                    xnddt = d2201 * std::cos(x2omi + xliL - g22) + d2211 * std::cos(xliL - g22)
                          + d3210 * std::cos(xomi + xliL - g32) + d3222 * std::cos(-xomi + xliL - g32)
                          + d5220 * std::cos(xomi + xliL - g52) + d5232 * std::cos(-xomi + xliL - g52)
                          + 2.0
                                * (d4410 * std::cos(x2omi + x2li - g44) + d4422 * std::cos(x2li - g44)
                                   + d5421 * std::cos(xomi + x2li - g54) + d5433 * std::cos(-xomi + x2li - g54));
                    xnddt *= xldot;
                }
                if (std::abs(t - atimeL) >= stepp) {
                    xliL   += xldot * delt + xndt * step2;
                    xniL   += xndt * delt + xnddt * step2;
                    atimeL += delt;
                } else {
                    ft = t - atimeL;
                    break;
                }
            }
            nm              = xniL + xndt * ft + xnddt * ft * ft * 0.5;
            const double xl = xliL + xldot * ft + xndt * ft * ft * 0.5;
            if (irez != 1) {
                mm = xl - 2.0 * nodem + 2.0 * theta;
            } else {
                mm = xl - nodem - argpm + theta;
            }
            dndt = nm - no_unkozai;
            nm   = no_unkozai + dndt;
        }
    }

    // ----------------------------------------------------------- initialization
    void init() {
        constexpr double temp4      = 1.5e-12;
        const double     ss         = 78.0 / radiusearthkm + 1.0;
        const double     qzms2ttemp = (120.0 - 78.0) / radiusearthkm;
        const double     qzms2t     = qzms2ttemp * qzms2ttemp * qzms2ttemp * qzms2ttemp;

        // initl
        const double eccsq  = ecco * ecco;
        const double omeosq = 1.0 - eccsq;
        const double rteosq = std::sqrt(omeosq);
        const double cosio  = std::cos(inclo);
        const double cosio2 = cosio * cosio;
        const double ak     = std::pow(xke / no_kozai, x2o3);
        const double d1     = 0.75 * j2 * (3.0 * cosio2 - 1.0) / (rteosq * omeosq);
        double       del    = d1 / (ak * ak);
        const double adel   = ak * (1.0 - del * del - del * (1.0 / 3.0 + 134.0 * del * del / 81.0));
        del                 = d1 / (adel * adel);
        no_unkozai          = no_kozai / (1.0 + del);
        const double ao     = std::pow(xke / no_unkozai, x2o3);
        const double sinio  = std::sin(inclo);
        const double po     = ao * omeosq;
        const double con42  = 1.0 - 5.0 * cosio2;
        con41               = -con42 - cosio2 - cosio2;
        const double posq   = po * po;
        const double rp     = ao * (1.0 - ecco);
        method              = 'n';
        if (opsmode == 'a') {
            const double     ts70  = epochDays - 7305.0;
            const double     ds70  = std::floor(ts70 + 1.0e-8);
            const double     tfrac = ts70 - ds70;
            constexpr double c1 = 1.72027916940703639e-2, thgr70 = 1.7321343856509374, fk5r = 5.07551419432269442e-15;
            const double     c1p2p = c1 + twopi;
            gsto                   = std::fmod(thgr70 + c1 * ds70 + c1p2p * tfrac + ts70 * ts70 * fk5r, twopi);
            if (gsto < 0.0) {
                gsto += twopi;
            }
        } else {
            gsto = gstime(epochDays + 2433281.5);
        }

        if (omeosq >= 0.0 || no_unkozai >= 0.0) {
            isimp = 0;
            if (rp < 220.0 / radiusearthkm + 1.0) {
                isimp = 1;
            }
            double       sfour  = ss;
            double       qzms24 = qzms2t;
            const double perige = (rp - 1.0) * radiusearthkm;
            if (perige < 156.0) {
                sfour = perige - 78.0;
                if (perige < 98.0) {
                    sfour = 20.0;
                }
                const double qzms24temp = (120.0 - sfour) / radiusearthkm;
                qzms24                  = qzms24temp * qzms24temp * qzms24temp * qzms24temp;
                sfour                   = sfour / radiusearthkm + 1.0;
            }
            const double pinvsq = 1.0 / posq;
            const double tsi    = 1.0 / (ao - sfour);
            eta                 = ao * ecco * tsi;
            const double etasq  = eta * eta;
            const double eeta   = ecco * eta;
            const double psisq  = std::abs(1.0 - etasq);
            const double coef   = qzms24 * std::pow(tsi, 4.0);
            const double coef1  = coef / std::pow(psisq, 3.5);
            const double cc2    = coef1 * no_unkozai
                             * (ao * (1.0 + 1.5 * etasq + eeta * (4.0 + etasq))
                                + 0.375 * j2 * tsi / psisq * con41 * (8.0 + 3.0 * etasq * (8.0 + etasq)));
            cc1        = bstar * cc2;
            double cc3 = 0.0;
            if (ecco > 1.0e-4) {
                cc3 = -2.0 * coef * tsi * j3oj2 * no_unkozai * sinio / ecco;
            }
            x1mth2 = 1.0 - cosio2;
            cc4    = 2.0 * no_unkozai * coef1 * ao * omeosq
                * (eta * (2.0 + 0.5 * etasq) + ecco * (0.5 + 2.0 * etasq)
                   - j2 * tsi / (ao * psisq)
                         * (-3.0 * con41 * (1.0 - 2.0 * eeta + etasq * (1.5 - 0.5 * eeta))
                            + 0.75 * x1mth2 * (2.0 * etasq - eeta * (1.0 + etasq)) * std::cos(2.0 * argpo)));
            cc5                 = 2.0 * coef1 * ao * omeosq * (1.0 + 2.75 * (etasq + eeta) + eeta * etasq);
            const double cosio4 = cosio2 * cosio2;
            const double temp1  = 1.5 * j2 * pinvsq * no_unkozai;
            const double temp2  = 0.5 * temp1 * j2 * pinvsq;
            const double temp3  = -0.46875 * j4 * pinvsq * pinvsq * no_unkozai;
            mdot                = no_unkozai + 0.5 * temp1 * rteosq * con41
                 + 0.0625 * temp2 * rteosq * (13.0 - 78.0 * cosio2 + 137.0 * cosio4);
            argpdot = -0.5 * temp1 * con42 + 0.0625 * temp2 * (7.0 - 114.0 * cosio2 + 395.0 * cosio4)
                    + temp3 * (3.0 - 36.0 * cosio2 + 49.0 * cosio4);
            const double xhdot1 = -temp1 * cosio;
            nodedot = xhdot1 + (0.5 * temp2 * (4.0 - 19.0 * cosio2) + 2.0 * temp3 * (3.0 - 7.0 * cosio2)) * cosio;
            const double xpidot = argpdot + nodedot;
            omgcof              = bstar * cc3 * std::cos(argpo);
            xmcof               = 0.0;
            if (ecco > 1.0e-4) {
                xmcof = -x2o3 * coef * bstar / eeta;
            }
            nodecf = 3.5 * omeosq * xhdot1 * cc1;
            t2cof  = 1.5 * cc1;
            if (std::abs(cosio + 1.0) > 1.5e-12) {
                xlcof = -0.25 * j3oj2 * sinio * (3.0 + 5.0 * cosio) / (1.0 + cosio);
            } else {
                xlcof = -0.25 * j3oj2 * sinio * (3.0 + 5.0 * cosio) / temp4;
            }
            aycof                  = -0.5 * j3oj2 * sinio;
            const double delmotemp = 1.0 + eta * std::cos(mo);
            delmo                  = delmotemp * delmotemp * delmotemp;
            sinmao                 = std::sin(mo);
            x7thm1                 = 7.0 * cosio2 - 1.0;

            // Deep space initialization.
            if (twopi / no_unkozai >= 225.0) {
                method             = 'd';
                isimp              = 1;
                const double tc    = 0.0;
                double       inclm = inclo;
                const DsCom  d     = dscom(epochDays, ecco, argpo, tc, inclo, nodeo, no_unkozai);
                e3                 = d.e3;
                ee2                = d.ee2;
                peo                = d.peo;
                pgho               = d.pgho;
                pho                = d.pho;
                pinco              = d.pinco;
                plo                = d.plo;
                se2                = d.se2;
                se3                = d.se3;
                sgh2               = d.sgh2;
                sgh3               = d.sgh3;
                sgh4               = d.sgh4;
                sh2                = d.sh2;
                sh3                = d.sh3;
                si2                = d.si2;
                si3                = d.si3;
                sl2                = d.sl2;
                sl3                = d.sl3;
                sl4                = d.sl4;
                xgh2               = d.xgh2;
                xgh3               = d.xgh3;
                xgh4               = d.xgh4;
                xh2                = d.xh2;
                xh3                = d.xh3;
                xi2                = d.xi2;
                xi3                = d.xi3;
                xl2                = d.xl2;
                xl3                = d.xl3;
                xl4                = d.xl4;
                zmol               = d.zmol;
                zmos               = d.zmos;
                dpper(0.0, true, ecco, inclo, nodeo, argpo, mo);
                double argpm = 0.0, nodem = 0.0, mm = 0.0, em = d.em, nm = d.nm, dndt = 0.0;
                dsinit(d, 0.0, tc, xpidot, em, argpm, inclm, mm, nm, nodem, dndt);
            }

            if (isimp != 1) {
                const double cc1sq = cc1 * cc1;
                d2                 = 4.0 * ao * tsi * cc1sq;
                const double temp  = d2 * tsi * cc1 / 3.0;
                d3                 = (17.0 * ao + sfour) * temp;
                d4                 = 0.5 * temp * ao * tsi * (221.0 * ao + 31.0 * sfour) * cc1;
                t3cof              = d2 + 2.0 * cc1sq;
                t4cof              = 0.25 * (3.0 * d3 + cc1 * (12.0 * d2 + 10.0 * cc1sq));
                t5cof = 0.2 * (3.0 * d4 + 12.0 * cc1 * d3 + 6.0 * d2 * d2 + 15.0 * cc1sq * (2.0 * d2 + cc1sq));
            }
        }
    }

    // ----------------------------------------------------------- propagation
    int sgp4(double t, StateVector& out) const {
        constexpr double temp4     = 1.5e-12;
        const double     vkmpersec = radiusearthkm * xke / 60.0;

        // Secular gravity and atmospheric drag.
        const double xmdf   = mo + mdot * t;
        const double argpdf = argpo + argpdot * t;
        const double nodedf = nodeo + nodedot * t;
        double       argpm  = argpdf;
        double       mm     = xmdf;
        const double t2     = t * t;
        double       nodem  = nodedf + nodecf * t2;
        double       tempa  = 1.0 - cc1 * t;
        double       tempe  = bstar * cc4 * t;
        double       templ  = t2cof * t2;
        if (isimp != 1) {
            const double delomg   = omgcof * t;
            const double delmtemp = 1.0 + eta * std::cos(xmdf);
            const double delm     = xmcof * (delmtemp * delmtemp * delmtemp - delmo);
            const double temp     = delomg + delm;
            mm                    = xmdf + temp;
            argpm                 = argpdf - temp;
            const double t3       = t2 * t;
            const double t4       = t3 * t;
            tempa                 = tempa - d2 * t2 - d3 * t3 - d4 * t4;
            tempe                 = tempe + bstar * cc5 * (std::sin(mm) - sinmao);
            templ                 = templ + t3cof * t3 + t4 * (t4cof + t * t5cof);
        }
        double nm    = no_unkozai;
        double em    = ecco;
        double inclm = inclo;
        if (method == 'd') {
            const double tc     = t;
            double       atimeL = atime, xliL = xli, xniL = xni, dndt = 0.0;
            dspace(t, tc, atimeL, em, argpm, inclm, xliL, mm, xniL, nodem, dndt, nm);
        }
        if (nm <= 0.0) {
            return 2;
        }
        const double am  = std::pow(xke / nm, x2o3) * tempa * tempa;
        nm               = xke / std::pow(am, 1.5);
        em              -= tempe;
        if (em >= 1.0 || em < -0.001) {
            return 1;
        }
        if (em < 1.0e-6) {
            em = 1.0e-6;
        }
        mm         += no_unkozai * templ;
        double xlm  = mm + argpm + nodem;
        nodem       = std::fmod(nodem, twopi);
        argpm       = std::fmod(argpm, twopi);
        xlm         = std::fmod(xlm, twopi);
        mm          = std::fmod(xlm - argpm - nodem, twopi);

        const double sinim = std::sin(inclm), cosim = std::cos(inclm);

        // Lunar-solar periodics.
        double ep = em, xincp = inclm, argpp = argpm, nodep = nodem, mp = mm;
        double sinip = sinim, cosip = cosim;
        double aycofL = aycof, xlcofL = xlcof;
        if (method == 'd') {
            dpper(t, false, ep, xincp, nodep, argpp, mp);
            if (xincp < 0.0) {
                xincp  = -xincp;
                nodep += pi;
                argpp -= pi;
            }
            if (ep < 0.0 || ep > 1.0) {
                return 3;
            }
            sinip  = std::sin(xincp);
            cosip  = std::cos(xincp);
            aycofL = -0.5 * j3oj2 * sinip;
            if (std::abs(cosip + 1.0) > 1.5e-12) {
                xlcofL = -0.25 * j3oj2 * sinip * (3.0 + 5.0 * cosip) / (1.0 + cosip);
            } else {
                xlcofL = -0.25 * j3oj2 * sinip * (3.0 + 5.0 * cosip) / temp4;
            }
        }

        // Long-period periodics.
        const double axnl = ep * std::cos(argpp);
        double       temp = 1.0 / (am * (1.0 - ep * ep));
        const double aynl = ep * std::sin(argpp) + temp * aycofL;
        const double xl   = mp + argpp + nodep + temp * xlcofL * axnl;

        // Solve Kepler's equation.
        const double u   = std::fmod(xl - nodep, twopi);
        double       eo1 = u, tem5 = 9999.9, sineo1 = 0.0, coseo1 = 0.0;
        for (int ktr = 1; std::abs(tem5) >= 1.0e-12 && ktr <= 10; ++ktr) {
            sineo1 = std::sin(eo1);
            coseo1 = std::cos(eo1);
            tem5   = 1.0 - coseo1 * axnl - sineo1 * aynl;
            tem5   = (u - aynl * coseo1 + axnl * sineo1 - eo1) / tem5;
            if (std::abs(tem5) >= 0.95) {
                tem5 = tem5 > 0.0 ? 0.95 : -0.95;
            }
            eo1 += tem5;
        }

        // Short-period preliminary quantities.
        const double ecose = axnl * coseo1 + aynl * sineo1;
        const double esine = axnl * sineo1 - aynl * coseo1;
        const double el2   = axnl * axnl + aynl * aynl;
        const double pl    = am * (1.0 - el2);
        if (pl < 0.0) {
            return 4;
        }
        const double rl     = am * (1.0 - ecose);
        const double rdotl  = std::sqrt(am) * esine / rl;
        const double rvdotl = std::sqrt(pl) / rl;
        const double betal  = std::sqrt(1.0 - el2);
        temp                = esine / (1.0 + betal);
        const double sinu   = am / rl * (sineo1 - aynl - axnl * temp);
        const double cosu   = am / rl * (coseo1 - axnl + aynl * temp);
        double       su     = std::atan2(sinu, cosu);
        const double sin2u  = (cosu + cosu) * sinu;
        const double cos2u  = 1.0 - 2.0 * sinu * sinu;
        temp                = 1.0 / pl;
        const double temp1  = 0.5 * j2 * temp;
        const double temp2  = temp1 * temp;

        double con41L = con41, x1mth2L = x1mth2, x7thm1L = x7thm1;
        if (method == 'd') {
            const double cosisq = cosip * cosip;
            con41L              = 3.0 * cosisq - 1.0;
            x1mth2L             = 1.0 - cosisq;
            x7thm1L             = 7.0 * cosisq - 1.0;
        }
        const double mrt    = rl * (1.0 - 1.5 * temp2 * betal * con41L) + 0.5 * temp1 * x1mth2L * cos2u;
        su                 -= 0.25 * temp2 * x7thm1L * sin2u;
        const double xnode  = nodep + 1.5 * temp2 * cosip * sin2u;
        const double xinc   = xincp + 1.5 * temp2 * cosip * sinip * cos2u;
        const double mvt    = rdotl - nm * temp1 * x1mth2L * sin2u / xke;
        const double rvdot  = rvdotl + nm * temp1 * (x1mth2L * cos2u + 1.5 * con41L) / xke;

        // Orientation vectors.
        const double sinsu = std::sin(su), cossu = std::cos(su);
        const double snod = std::sin(xnode), cnod = std::cos(xnode);
        const double sini = std::sin(xinc), cosi = std::cos(xinc);
        const double xmx = -snod * cosi;
        const double xmy = cnod * cosi;
        const Vec3   uv{xmx * sinsu + cnod * cossu, xmy * sinsu + snod * cossu, sini * sinsu};
        const Vec3   vv{xmx * cossu - cnod * sinsu, xmy * cossu - snod * sinsu, sini * cossu};
        out.r = uv * (mrt * radiusearthkm * 1000.0);
        out.v = (uv * mvt + vv * rvdot) * (vkmpersec * 1000.0);
        if (mrt < 1.0) {
            return 6;
        }
        return 0;
    }
};

Sgp4::Sgp4(const Tle& tle, Sgp4Gravity gravity, Sgp4OpsMode mode)
    : impl_(std::make_unique<Impl>()) {
    // Check what the theory uses (the two derivatives of the mean motion are carried but not used).
    requireFinite(tle.epoch.jd().jd1 + tle.epoch.jd().jd2, "epoch");
    requireFinite(tle.bstar, "B*");
    requireFinite(tle.inclination, "inclination");
    requireFinite(tle.raan, "right ascension of the ascending node");
    requireFinite(tle.eccentricity, "eccentricity");
    requireFinite(tle.argPerigee, "argument of perigee");
    requireFinite(tle.meanAnomaly, "mean anomaly");
    requireFinite(tle.meanMotion, "mean motion");
    // The reference code tests for error 2 after the Kozai correction, which takes a fractional power of the mean
    // motion: a negative one has turned into NaN by then and is not caught.
    if (tle.meanMotion < 0.0) {
        throw Sgp4Error(2, std::string(errorMessage(2)) + " (at epoch)");
    }
    Impl& s     = *impl_;
    s.tle       = tle;
    s.epochTime = tle.epoch;
    s.opsmode   = mode == Sgp4OpsMode::Afspc ? 'a' : 'i';
    s.setGravity(gravity);
    s.bstar             = tle.bstar;
    s.ndot              = tle.meanMotionDot / (xpdotp * 1440.0);
    s.nddot             = tle.meanMotionDdot / (xpdotp * 1440.0 * 1440.0);
    s.ecco              = tle.eccentricity;
    s.argpo             = tle.argPerigee;
    s.inclo             = tle.inclination;
    s.mo                = tle.meanAnomaly;
    s.no_kozai          = tle.meanMotion / xpdotp;
    s.nodeo             = tle.raan;
    const JulianDate jd = tle.epoch.to(TimeScale::UTC).jd();
    // The reference implementation counts the epoch in days of 86400 s, straight from the day fraction of the
    // element set, so the fraction of a UTC day that ends with a leap second is put back on that footing (the
    // factor is exactly 1 on any other day).
    const double dayFraction = jd.jd2 / tleDayPerUtcDay(jd.jd1);
    // Deliberately reproduce the reference implementation's arithmetic, which sums the two JD parts into one
    // double before subtracting (jdsatepoch + jdsatepochF - 2433281.5). The resulting ~1e-10 day rounding
    // reaches the lunar-solar terms, and matching it gives sub-millimeter agreement with Vallado's test output
    // even for extreme deep-space orbits.
    s.epochDays = (jd.jd1 + dayFraction) - 2433281.5;
    s.init();
    // Check the epoch state as Vallado does (errors 1-4 at epoch are fatal).
    StateVector tmp;
    const int   code = s.sgp4(0.0, tmp);
    if (code != 0 && code != 6) {
        throw Sgp4Error(code, std::string(errorMessage(code)) + " (at epoch)");
    }
    // Finite elements can still overflow on the way (a drag term of 1e300, say). No range test catches that,
    // because the tests are false for NaN, so look at the state itself.
    for (const double component : {tmp.r.x, tmp.r.y, tmp.r.z, tmp.v.x, tmp.v.y, tmp.v.z}) {
        if (!std::isfinite(component)) {
            throw Error("Sgp4: the elements give a non-finite state at the epoch");
        }
    }
}

Sgp4::~Sgp4() = default;
Sgp4::Sgp4(const Sgp4& o)
    : impl_(std::make_unique<Impl>(*o.impl_)) {}
Sgp4& Sgp4::operator=(const Sgp4& o) {
    if (this != &o) {
        impl_ = std::make_unique<Impl>(*o.impl_);
    }
    return *this;
}
Sgp4::Sgp4(Sgp4&&) noexcept            = default;
Sgp4& Sgp4::operator=(Sgp4&&) noexcept = default;

int Sgp4::propagate(double minutes, StateVector& out) const noexcept {
    return impl_->sgp4(minutes, out);
}

StateVector Sgp4::propagate(double minutes) const {
    StateVector out;
    const int   code = impl_->sgp4(minutes, out);
    if (code != 0) {
        throw Sgp4Error(code, errorMessage(code));
    }
    return out;
}

StateVector Sgp4::propagate(const Time& t) const {
    return propagate((t - impl_->epochTime) / 60.0);
}

const Tle& Sgp4::tle() const noexcept {
    return impl_->tle;
}
const Time& Sgp4::epoch() const noexcept {
    return impl_->epochTime;
}
bool Sgp4::isDeepSpace() const noexcept {
    return impl_->method == 'd';
}

std::ostream& operator<<(std::ostream& os, Sgp4Gravity value) {
    return detail::writeEnum(os, "Sgp4Gravity", value);
}

std::ostream& operator<<(std::ostream& os, Sgp4OpsMode value) {
    return detail::writeEnum(os, "Sgp4OpsMode", value);
}

} // namespace astro
