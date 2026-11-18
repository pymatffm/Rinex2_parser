"""Synthetic RINEX 3.04 / 4.00 generator (test fixtures; NOT real receiver data)."""
import numpy as np
from datetime import datetime, timedelta

XYZ = (4075539.8920, 931735.2680, 4801629.4330)   # Wettzell-like
OBS = {"G": ["C1C", "L1C", "D1C", "S1C", "C2W", "L2W", "S2W", "C5Q", "L5Q", "S5Q"],
       "E": ["C1C", "L1C", "S1C", "C5Q", "L5Q", "S5Q", "C7Q", "L7Q", "S7Q"],
       "C": ["C2I", "L2I", "S2I", "C7I", "L7I", "S7I", "C6I", "L6I", "S6I"]}
WL = {"1": .1902937, "2": .2442103, "5": .2548280, "7": .2483493, "6": .2365}
COUNTS = {"G": 32, "E": 36, "C": 46}

def hdr(line, label): return f"{line:<60}{label}\n"

def obs_header(ver, tstart):
    s = hdr(f"{ver:9.2f}           OBSERVATION DATA    M: MIXED", "RINEX VERSION / TYPE")
    s += hdr("rinex_tool synthetic", "PGM / RUN BY / DATE")
    s += hdr("SYNT", "MARKER NAME")
    s += hdr("SYNT-RX             SIMULATOR           1.0", "REC # / TYPE / VERS")
    s += hdr("0                   TRM59800.00     NONE", "ANT # / TYPE")
    s += hdr("%14.4f%14.4f%14.4f" % XYZ, "APPROX POSITION XYZ")
    s += hdr("        0.0000        0.0000        0.0000", "ANTENNA: DELTA H/E/N")
    for sy, codes in OBS.items():   # 13 codes per line
        for k in range(0, len(codes), 13):
            chunk = codes[k:k+13]
            head = f"{sy}  {len(codes):3d}" if k == 0 else " " * 6
            s += hdr(head + "".join(f" {c}" for c in chunk), "SYS / # / OBS TYPES")
    s += hdr("DBHZ", "SIGNAL STRENGTH UNIT")
    s += hdr("    30.000", "INTERVAL")
    t = tstart
    s += hdr(f"  {t.year:4d}    {t.month:2d}    {t.day:2d}    {t.hour:2d}    {t.minute:2d}   {t.second:2d}.0000000     GPS", "TIME OF FIRST OBS")
    s += hdr("    18", "LEAP SECONDS")
    s += hdr("", "END OF HEADER")
    return s

def make_obs(path, ver=3.04, hours=2.0, seed=1):
    rng = np.random.default_rng(seed)
    t0 = datetime(2024, 3, 14, 0, 0, 0)
    n = int(hours * 120)
    sats = {}
    for sy, c in COUNTS.items():
        for p in range(1, c + 1):
            ph = rng.uniform(0, 2 * np.pi); per = rng.uniform(40, 100)
            sats[f"{sy}{p:02d}"] = (ph, per, rng.uniform(.3, .6))
    out = [obs_header(ver, t0)]
    for k in range(n):
        t = t0 + timedelta(seconds=30 * k)
        rows = []
        for sat, (ph, per, thr) in sats.items():
            elev = np.sin(2 * np.pi * k / (per * 2) + ph)
            if elev < thr - .3: continue
            rows.append((sat, elev))
        out.append(f"> {t.year:4d} {t.month:02d} {t.day:02d} {t.hour:02d} {t.minute:02d} {t.second:2d}.0000000  0{len(rows):3d}\n")
        for sat, elev in rows:
            sy = sat[0]; line = sat
            rng_m = 2.2e7 + 4e6 * (1 - elev) + 30 * k
            for code in OBS[sy]:
                ty, b = code[0], code[1]
                if ty == "C": v = rng_m + rng.normal(0, .3)
                elif ty == "L": v = rng_m / WL[b] + rng.normal(0, .01)
                elif ty == "D": v = -2000 + 3 * k
                else: v = 32 + 16 * max(elev, 0) + (b == "5") * 2 + rng.normal(0, 1)
                lli = "1" if (ty == "L" and rng.random() < .002) else " "
                line += f"{v:14.3f}{lli} " if ty != "S" else f"{v:14.3f}  "
            out.append(line.rstrip() + "\n")
    open(path, "w").write("".join(out))

def _e(v): return f"{v:19.12e}"
def eph_lines(sat, t, rng, ver):
    l = [f"{sat} {t.year:4d} {t.month:02d} {t.day:02d} {t.hour:02d} {t.minute:02d} {t.second:02d}"
         + _e(1.3e-4 * rng.normal()) + _e(7e-12) + _e(0.0)]
    vals = [[98, -1.7, 4.6e-9, 2.1], [-1.8e-7, 3.3e-4, 8.2e-6, 5153.8],
            [t.hour * 3600 + 7200, 5e-8, -1.2, 9e-8], [.96, 200.5, -1.5, -8e-9],
            [-1e-10, 1, 2300, 0], [2.0, 0, 5e-9, 98], [t.hour * 3600 + 100, 4.0]]
    for v in vals: l.append("    " + "".join(_e(x) for x in v))
    return l

def make_nav(path, ver=3.04, seed=2):
    rng = np.random.default_rng(seed)
    t = datetime(2024, 3, 14, 2, 0, 0)
    s = hdr(f"{ver:9.2f}           N: GNSS NAV DATA    M: MIXED", "RINEX VERSION / TYPE")
    if ver < 4:
        s += hdr("GPSA   1.1176E-08  2.2352E-08 -5.9605E-08 -1.1921E-07", "IONOSPHERIC CORR")
        s += hdr("GPUT  2.synthetic".replace("2.synthetic", " 1.8626451492E-09 1.776357E-15 589824 2305"), "TIME SYSTEM CORR")
    s += hdr("    18", "LEAP SECONDS"); s += hdr("", "END OF HEADER")
    for sat, msg in (("G04", "LNAV"), ("E11", "INAV"), ("C06", "D1D2")):
        if ver >= 4: s += f"> EPH {sat} {msg}\n"
        s += "\n".join(eph_lines(sat, t, rng, ver)) + "\n"
    if ver >= 4:
        s += "> STO E   IFNV\n    2024 03 14 00 00 00 GAUT                                  UTCGAL\n"
        s += "     1.735000000000e+05-1.862645149231e-09 0.000000000000e+00 0.000000000000e+00 0.000000000000e+00\n"
        s += "> ION G14 CNVX\n    2024 03 14 23 30 42 7.450580596924e-09 2.235174179077e-08-5.960464477539e-08\n"
        s += "    -1.192092895508e-07 9.216000000000e+04 1.310720000000e+05-6.553600000000e+04\n"
        s += "    -5.242880000000e+05\n"
    open(path, "w").write(s)
