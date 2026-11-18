#!/usr/bin/env python3
"""
rinex_tool - RINEX 3.x / 4.x observation and navigation parser, time-scale
conversion, tidy pandas extraction, summary statistics and plots.

Supported: GPS (G), Galileo (E), BeiDou (C) as first-class systems; R, J, I, S are
parsed generically. RINEX 2 is rejected with a clear error.

Time scales (all continuous except UTC):
    GPST = TAI - 19 s          GST  = GPST (aligned, different week rollover)
    BDT  = GPST - 14 s         UTC  = TAI - (TAI-UTC) from leap-second table
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

SYSTEMS = {"G": "GPS", "E": "Galileo", "C": "BeiDou", "R": "GLONASS",
           "J": "QZSS", "I": "NavIC", "S": "SBAS"}

# ---------------------------------------------------------------- time scales
# (UTC date the offset takes effect, TAI-UTC in seconds after that instant)
_LEAP = [("1980-01-01", 19), ("1981-07-01", 20), ("1982-07-01", 21),
         ("1983-07-01", 22), ("1985-07-01", 23), ("1988-01-01", 24),
         ("1990-01-01", 25), ("1991-01-01", 26), ("1992-07-01", 27),
         ("1993-07-01", 28), ("1994-07-01", 29), ("1996-01-01", 30),
         ("1997-07-01", 31), ("1999-01-01", 32), ("2006-01-01", 33),
         ("2009-01-01", 34), ("2012-07-01", 35), ("2015-07-01", 36),
         ("2017-01-01", 37)]
_LEAP_UTC = np.array([np.datetime64(d, "ns") for d, _ in _LEAP])
_LEAP_VAL = np.array([v for _, v in _LEAP])
# the same instants expressed on the TAI scale
_LEAP_TAI = _LEAP_UTC + (_LEAP_VAL * 1_000_000_000).astype("timedelta64[ns]")

GPS_EPOCH = pd.Timestamp("1980-01-06")
_S = lambda x: pd.to_timedelta(x, unit="s")  # noqa: E731

#: scale name -> offset such that  scale = GPST + offset  (fixed-offset scales)
_FIXED = {"GPS": 0, "GAL": 0, "GST": 0, "GPST": 0, "QZS": 0, "IRN": 0,
          "BDT": -14, "BDS": -14, "TAI": 19}


def tai_minus_utc_from_utc(utc: pd.DatetimeIndex) -> np.ndarray:
    i = np.searchsorted(_LEAP_UTC, utc.values, side="right") - 1
    if (i < 0).any():
        raise ValueError("UTC time before 1980-01-01 not supported")
    return _LEAP_VAL[i]


def tai_minus_utc_from_tai(tai: pd.DatetimeIndex) -> np.ndarray:
    i = np.searchsorted(_LEAP_TAI, tai.values, side="right") - 1
    if (i < 0).any():
        raise ValueError("TAI time before 1980-01-01 not supported")
    return _LEAP_VAL[i]


def _norm(scale: str) -> str:
    s = scale.upper().strip()
    return "UTC" if s in ("GLO", "UTC") else s


def to_gpst(t, scale: str) -> pd.DatetimeIndex:
    """Convert timestamps expressed in `scale` to GPST."""
    t = pd.DatetimeIndex(t)
    s = _norm(scale)
    if s == "UTC":
        return t + pd.to_timedelta(tai_minus_utc_from_utc(t) - 19, unit="s")
    if s in _FIXED:
        return t - _S(_FIXED[s])
    raise ValueError(f"unknown time scale {scale!r}")


def from_gpst(t, scale: str) -> pd.DatetimeIndex:
    """Convert GPST timestamps to `scale`."""
    t = pd.DatetimeIndex(t)
    s = _norm(scale)
    if s == "UTC":
        tai = t + _S(19)
        return tai - pd.to_timedelta(tai_minus_utc_from_tai(tai), unit="s")
    if s in _FIXED:
        return t + _S(_FIXED[s])
    raise ValueError(f"unknown time scale {scale!r}")


def convert_time(t, src: str, dst: str) -> pd.DatetimeIndex:
    """Convert between GPST / GST / BDT / UTC / TAI (any pair)."""
    return from_gpst(to_gpst(t, src), dst)


def gps_week_tow(gpst) -> pd.DataFrame:
    """GPS week number (continuous, no rollover) and time-of-week for GPST times."""
    d = pd.DatetimeIndex(gpst) - GPS_EPOCH
    sec = d.total_seconds().to_numpy()
    wk = np.floor(sec / 604800).astype(int)
    return pd.DataFrame({"week": wk, "tow": sec - wk * 604800})


# ------------------------------------------------------------------- helpers
def _label(line: str) -> str:
    return line[60:].strip()


def _f(s: str) -> float:
    s = s.strip().replace("D", "E").replace("d", "e")
    return float(s) if s else np.nan


def _split_header_time(y, mo, d, h, mi, sec) -> pd.Timestamp:
    return (pd.Timestamp(int(y), int(mo), int(d), int(h), int(mi))
            + pd.Timedelta(seconds=float(sec)))


def _check_version(line: str) -> tuple[float, str, str]:
    ver = float(line[0:9])
    if ver < 3.0:
        raise ValueError(f"RINEX {ver} is not supported (need 3.x or 4.x)")
    ftype = line[20:21]
    sat_sys = line[40:41].strip() or "G"
    return ver, ftype, sat_sys


# ------------------------------------------------------------- observations
@dataclass
class ObsHeader:
    version: float = 0.0
    file_type: str = "O"
    marker_name: str = ""
    receiver: str = ""
    antenna: str = ""
    approx_xyz: Optional[tuple] = None
    antenna_delta_hen: Optional[tuple] = None
    obs_types: Dict[str, List[str]] = field(default_factory=dict)
    time_first_obs: Optional[pd.Timestamp] = None
    time_last_obs: Optional[pd.Timestamp] = None
    time_system: str = "GPS"
    interval: Optional[float] = None
    leap_seconds: Optional[int] = None
    signal_strength_unit: str = ""
    comments: List[str] = field(default_factory=list)

    @property
    def approx_llh(self):
        """Geodetic (lat deg, lon deg, h m) on WGS-84 from the header XYZ."""
        if self.approx_xyz is None:
            return None
        return xyz_to_llh(*self.approx_xyz)


def xyz_to_llh(x: float, y: float, z: float) -> tuple:
    a, f = 6378137.0, 1 / 298.257223563
    e2 = f * (2 - f)
    lon = np.arctan2(y, x)
    p = np.hypot(x, y)
    if p < 1e-9:  # pole / origin
        return (90.0 if z >= 0 else -90.0, 0.0, abs(z) - a * (1 - f))
    lat = np.arctan2(z, p * (1 - e2))
    for _ in range(10):
        n = a / np.sqrt(1 - e2 * np.sin(lat) ** 2)
        h = p / np.cos(lat) - n
        lat = np.arctan2(z, p * (1 - e2 * n / (n + h)))
    n = a / np.sqrt(1 - e2 * np.sin(lat) ** 2)
    h = p / np.cos(lat) - n
    return float(np.degrees(lat)), float(np.degrees(lon)), float(h)


def parse_obs_header(lines: List[str]) -> tuple[ObsHeader, int]:
    h = ObsHeader()
    cur_sys = None
    expect = 0
    for i, line in enumerate(lines):
        lab = _label(line)
        if lab == "RINEX VERSION / TYPE":
            h.version, h.file_type, _ = _check_version(line)
            if h.file_type != "O":
                raise ValueError(f"not an observation file (type {h.file_type!r})")
        elif lab == "MARKER NAME":
            h.marker_name = line[:60].strip()
        elif lab == "REC # / TYPE / VERS":
            h.receiver = line[20:40].strip()
        elif lab == "ANT # / TYPE":
            h.antenna = line[20:40].strip()
        elif lab == "APPROX POSITION XYZ":
            h.approx_xyz = tuple(float(v) for v in line[:60].split()[:3])
        elif lab == "ANTENNA: DELTA H/E/N":
            h.antenna_delta_hen = tuple(float(v) for v in line[:60].split()[:3])
        elif lab == "SYS / # / OBS TYPES":
            if line[0] != " ":  # new constellation block
                cur_sys = line[0]
                expect = int(line[3:6])
                h.obs_types[cur_sys] = []
            codes = line[7:60].split()
            h.obs_types[cur_sys].extend(codes)
        elif lab == "TIME OF FIRST OBS":
            t = line[:43].split()
            h.time_first_obs = _split_header_time(*t[:6])
            sysname = line[48:51].strip()
            h.time_system = sysname or {"G": "GPS", "E": "GAL", "C": "BDT",
                                        "R": "GLO"}.get(_sat_sys_of(h), "GPS")
        elif lab == "TIME OF LAST OBS":
            t = line[:43].split()
            h.time_last_obs = _split_header_time(*t[:6])
        elif lab == "INTERVAL":
            h.interval = float(line[:10])
        elif lab == "LEAP SECONDS":
            h.leap_seconds = int(line[:6])
        elif lab == "SIGNAL STRENGTH UNIT":
            h.signal_strength_unit = line[:20].strip()
        elif lab == "COMMENT":
            h.comments.append(line[:60].strip())
        elif lab == "END OF HEADER":
            return h, i + 1
    raise ValueError("END OF HEADER not found")


def _sat_sys_of(h: ObsHeader) -> str:
    return next(iter(h.obs_types), "G")


_OBS_RE = re.compile(r"^([CLDSIX])(\d)(\w?)$")


def _split_code(code: str) -> tuple[str, str, str]:
    m = _OBS_RE.match(code)
    return (m.group(1), m.group(2), m.group(3)) if m else (code[:1], code[1:2], code[2:])


def read_obs(path, systems: str = "GEC"):
    """Parse a RINEX 3/4 observation file.

    Returns (header, df). `df` is tidy/long: one row per satellite, epoch, signal.

    Columns: time (file scale), time_gpst, time_utc, time_tai, sat, system, prn,
    code (e.g. C1C), type (C pseudorange, L carrier phase [cycles], D Doppler,
    S signal strength [dB-Hz]), band, attr, value, lli, ssi.
    """
    with open(path, "r", errors="replace") as fh:
        lines = fh.read().splitlines()
    header, i = parse_obs_header(lines)
    keep = set(systems)
    times, recs = [], []
    n = len(lines)
    while i < n:
        line = lines[i]
        if not line.startswith(">"):
            i += 1
            continue
        tok = line[1:].split()
        t = _split_header_time(*tok[:6])
        flag, nsat = int(tok[6]), int(tok[7])
        i += 1
        if flag > 1:                    # event record: skip text lines
            i += nsat
            continue
        for _ in range(nsat):
            row = lines[i]
            i += 1
            s = row[:3].replace(" ", "0")
            sys_ = s[0]
            if sys_ not in keep or sys_ not in header.obs_types:
                continue
            codes = header.obs_types[sys_]
            row = row[3:].ljust(16 * len(codes))
            for k, code in enumerate(codes):
                fld = row[16 * k:16 * k + 16]
                v = fld[:14].strip()
                if not v:
                    continue
                lli, ssi = fld[14:15].strip(), fld[15:16].strip()
                recs.append((t, s, sys_, int(s[1:]), code, float(v),
                             int(lli) if lli else 0, int(ssi) if ssi else 0))
    cols = ["time", "sat", "system", "prn", "code", "value", "lli", "ssi"]
    df = pd.DataFrame(recs, columns=cols)
    if df.empty:
        return header, df
    parts = df["code"].map(_split_code)
    df["type"] = [p[0] for p in parts]
    df["band"] = [p[1] for p in parts]
    df["attr"] = [p[2] for p in parts]
    df["time"] = pd.to_datetime(df["time"])
    uniq = pd.DatetimeIndex(df["time"].unique())
    gp = to_gpst(uniq, header.time_system)
    maps = {"time_gpst": gp, "time_utc": from_gpst(gp, "UTC"),
            "time_tai": from_gpst(gp, "TAI")}
    for name, idx in maps.items():
        df[name] = df["time"].map(dict(zip(uniq, idx)))
    df = df[["time", "time_gpst", "time_utc", "time_tai", "sat", "system", "prn",
             "code", "type", "band", "attr", "value", "lli", "ssi"]]
    return header, df


def observations_wide(df: pd.DataFrame, types=("C", "L")) -> pd.DataFrame:
    """Pivot to one row per (time, sat) with a column per observation code."""
    sub = df[df["type"].isin(types)]
    w = sub.pivot_table(index=["time", "time_gpst", "time_utc", "time_tai", "sat"],
                        columns="code", values="value", aggfunc="first")
    return w.reset_index()


def carrier_phase_and_pseudorange(df: pd.DataFrame) -> pd.DataFrame:
    """Only the C (pseudorange, m) and L (carrier phase, cycles) rows."""
    return df[df["type"].isin(["C", "L"])].reset_index(drop=True)


# ---------------------------------------------------------------- navigation
_NAMES_COMMON = [
    ["iode", "crs", "delta_n", "m0"],
    ["cuc", "e", "cus", "sqrt_a"],
    ["toe", "cic", "omega0", "cis"],
    ["i0", "crc", "omega", "omega_dot"],
    ["idot", "codes_l2", "week", "l2p_flag"],
    ["sv_accuracy", "sv_health", "tgd_bgd1", "iodc_bgd2"],
    ["t_transmission", "fit_interval"],
]
_NAMED_MSG = {"LNAV", "INAV", "FNAV", "D1", "D2", "D1D2", "FDMA"}


@dataclass
class NavData:
    version: float
    header: dict
    ephemeris: pd.DataFrame
    sto: pd.DataFrame
    ion: pd.DataFrame


def _fields(line: str, start: int) -> List[float]:
    out = []
    for k in range(start, len(line), 19):
        seg = line[k:k + 19]
        if seg.strip():
            out.append(_f(seg))
    return out


def _parse_eph(rec: List[str], msg: Optional[str], toc_scale: str) -> dict:
    first = rec[0]
    sat = first[:3].replace(" ", "0")
    t = [int(x) for x in first[4:23].split()]
    toc = pd.Timestamp(*t)
    vals = _fields(first, 23)
    row = {"sat": sat, "system": sat[0], "msg_type": msg or "LNAV",
           "toc": toc, "af0": vals[0], "af1": vals[1], "af2": vals[2]}
    flat = []
    for ln in rec[1:]:
        flat.extend(_fields(ln, 4))
    row["raw"] = flat
    if sat[0] in "GECJI" and (msg is None or msg in _NAMED_MSG):
        k = 0
        for names in _NAMES_COMMON:
            for nm in names:
                if k < len(flat):
                    row[nm] = flat[k]
                k += 1
    row["toc_scale"] = toc_scale
    return row


_TIME_SCALE_OF_SYS = {"G": "GPS", "E": "GAL", "C": "BDT", "J": "GPS",
                      "I": "GPS", "R": "UTC", "S": "GPS"}


def read_nav(path, systems: str = "GEC") -> NavData:
    """Parse a RINEX 3/4 navigation file (v4 EPH/STO/ION records supported)."""
    with open(path, "r", errors="replace") as fh:
        lines = fh.read().splitlines()
    hdr: dict = {"iono": {}, "time_corr": {}}
    i = 0
    for i, line in enumerate(lines):
        lab = _label(line)
        if lab == "RINEX VERSION / TYPE":
            ver, ft, _ = _check_version(line)
            if ft != "N":
                raise ValueError(f"not a navigation file (type {ft!r})")
            hdr["version"] = ver
        elif lab == "IONOSPHERIC CORR":
            tok = line[:60].split()
            hdr["iono"][tok[0]] = [_f(v) for v in tok[1:5]]
        elif lab == "TIME SYSTEM CORR":
            tok = line[:60].split()
            if len(tok) >= 4:
                hdr["time_corr"][tok[0]] = {"a0": _f(tok[1]), "a1": _f(tok[2]),
                                            "t": int(float(tok[3])),
                                            "w": int(float(tok[4])) if len(tok) > 4 else 0}
        elif lab == "LEAP SECONDS":
            hdr["leap_seconds"] = int(line[:6])
        elif lab == "END OF HEADER":
            i += 1
            break
    else:
        raise ValueError("END OF HEADER not found")
    ver = hdr["version"]
    eph, sto, ion = [], [], []
    keep = set(systems)
    body = lines[i:]

    # group into records
    records: List[tuple] = []        # (kind, msg, sat_hint, [lines])
    cur = None
    for ln in body:
        if not ln.strip():
            continue
        if ln.startswith(">"):
            tk = ln[1:].split()
            cur = [tk[0], tk[-1], tk[1] if len(tk) > 2 else "", []]
            records.append(cur)
        elif ver < 4 and ln[0] != " ":     # v3: new satellite record
            cur = ["EPH", None, ln[:3], [ln]]
            records.append(cur)
        elif cur is not None:
            cur[3].append(ln)
    for kind, msg, sat_hint, rec in records:
        if not rec:
            continue
        if kind == "EPH":
            sat = rec[0][:3].replace(" ", "0")
            if sat[0] not in keep:
                continue
            eph.append(_parse_eph(rec, msg, _TIME_SCALE_OF_SYS.get(sat[0], "GPS")))
        elif kind == "STO":
            e = rec[0]
            nums = _fields(rec[1], 4) if len(rec) > 1 else []
            sto.append({"system": sat_hint[:1], "msg_type": msg,
                        "epoch": pd.Timestamp(*[int(x) for x in e[4:23].split()]),
                        "label": e[24:28], "t_tm": nums[0] if nums else np.nan,
                        "a0": nums[1] if len(nums) > 1 else np.nan,
                        "a1": nums[2] if len(nums) > 2 else np.nan,
                        "a2": nums[3] if len(nums) > 3 else np.nan})
        elif kind == "ION":
            e = rec[0]
            nums = _fields(e, 23)
            for ln in rec[1:]:
                nums.extend(_fields(ln, 4))
            ion.append({"sat": sat_hint, "system": sat_hint[:1], "msg_type": msg,
                        "epoch": pd.Timestamp(*[int(x) for x in e[4:23].split()]),
                        "params": nums})
    e_df = pd.DataFrame(eph)
    if not e_df.empty:
        # clock epochs -> GPST for all scales
        e_df["toc_gpst"] = pd.NaT
        for sc, idx in e_df.groupby("toc_scale").groups.items():
            e_df.loc[idx, "toc_gpst"] = to_gpst(e_df.loc[idx, "toc"], sc)
    return NavData(ver, hdr, e_df, pd.DataFrame(sto), pd.DataFrame(ion))


# ------------------------------------------------------------ stats & plots
_BAND_LABEL = {
    "G": {"1": "L1", "2": "L2", "5": "L5"},
    "E": {"1": "E1", "5": "E5a", "7": "E5b", "8": "E5", "6": "E6"},
    "C": {"1": "B1C", "2": "B1I", "5": "B2a", "7": "B2b", "8": "B2", "6": "B3"},
}


def band_label(system: str, band: str) -> str:
    return _BAND_LABEL.get(system, {}).get(band, f"{system}{band}")


def summarize(header: ObsHeader, df: pd.DataFrame) -> dict:
    out: dict = {
        "rinex_version": header.version, "marker": header.marker_name,
        "receiver": header.receiver, "antenna": header.antenna,
        "approx_xyz_m": header.approx_xyz,
        "approx_llh": header.approx_llh, "time_system": header.time_system,
        "obs_types": header.obs_types,
    }
    if df.empty:
        return out
    ep = np.sort(df["time"].unique())
    dt = np.diff(ep).astype("timedelta64[ms]").astype(float) / 1000
    out.update({
        "first_epoch_gpst": str(df["time_gpst"].min()),
        "last_epoch_gpst": str(df["time_gpst"].max()),
        "first_epoch_utc": str(df["time_utc"].min()),
        "first_epoch_tai": str(df["time_tai"].min()),
        "n_epochs": int(len(ep)),
        "median_interval_s": float(np.median(dt)) if len(dt) else None,
        "n_gaps": int((dt > 1.5 * np.median(dt)).sum()) if len(dt) else 0,
        "n_rows": int(len(df)),
    })
    per_sys = {}
    for s, g in df.groupby("system"):
        per_sys[s] = {"name": SYSTEMS.get(s, s), "n_satellites": int(g["sat"].nunique()),
                      "satellites": sorted(g["sat"].unique().tolist()),
                      "n_pseudorange": int((g["type"] == "C").sum()),
                      "n_carrier_phase": int((g["type"] == "L").sum()),
                      "cycle_slip_flags": int(((g["type"] == "L") & (g["lli"] & 1 > 0)).sum())}
    out["per_system"] = per_sys
    return out


def snr_table(df: pd.DataFrame) -> pd.DataFrame:
    s = df[df["type"] == "S"]
    if s.empty:
        return pd.DataFrame()
    g = s.groupby(["system", "code"])["value"].agg(
        n="count", mean="mean", std="std", min="min", median="median", max="max")
    g = g.reset_index()
    g.insert(2, "signal", [band_label(a, c[1]) for a, c in zip(g["system"], g["code"])])
    return g.round(2)


def visibility(df: pd.DataFrame) -> pd.DataFrame:
    """Satellites tracked per epoch and system (any observation present)."""
    v = df.drop_duplicates(["time", "sat"])
    t = v.groupby(["time", "system"])["sat"].nunique().unstack(fill_value=0)
    t["total"] = t.sum(axis=1)
    return t


def make_plots(df: pd.DataFrame, outdir: Path, title: str = "") -> List[Path]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates

    outdir.mkdir(parents=True, exist_ok=True)
    colors = {"G": "#1f77b4", "E": "#2ca02c", "C": "#d62728", "R": "#9467bd",
              "J": "#8c564b", "I": "#e377c2", "S": "#7f7f7f"}
    files = []
    fmt = mdates.DateFormatter("%H:%M")

    # 1. SNR vs time, mean per system/band
    s = df[df["type"] == "S"].copy()
    if not s.empty:
        s["signal"] = [band_label(a, b) for a, b in zip(s["system"], s["band"])]
        fig, ax = plt.subplots(figsize=(11, 5))
        for (sy, sig), g in s.groupby(["system", "signal"]):
            m = g.groupby("time_gpst")["value"].mean()
            ax.plot(m.index, m.values, lw=1.2, label=f"{SYSTEMS.get(sy, sy)} {sig}",
                    color=colors.get(sy), ls=["-", "--", ":", "-."][hash(sig) % 4])
        ax.set(xlabel="GPST", ylabel="Mean C/N0 (dB-Hz)",
               title=f"Signal-to-noise ratio over time {title}".strip())
        ax.xaxis.set_major_formatter(fmt)
        ax.grid(alpha=.3)
        ax.legend(ncol=3, fontsize=8)
        p = outdir / "snr_vs_time.png"
        fig.tight_layout()
        fig.savefig(p, dpi=150)
        plt.close(fig)
        files.append(p)

        # 2. SNR distribution per system/band
        fig, ax = plt.subplots(figsize=(11, 5))
        groups = [(f"{sy}\n{sig}", g["value"].values)
                  for (sy, sig), g in s.groupby(["system", "signal"])]
        ax.boxplot([g for _, g in groups], tick_labels=[n for n, _ in groups],
                   showfliers=False)
        ax.set(ylabel="C/N0 (dB-Hz)", title="SNR distribution by signal")
        ax.grid(alpha=.3, axis="y")
        p = outdir / "snr_distribution.png"
        fig.tight_layout()
        fig.savefig(p, dpi=150)
        plt.close(fig)
        files.append(p)

    # 3. satellites visible vs time
    vis = visibility(df)
    fig, ax = plt.subplots(figsize=(11, 5))
    syss = [c for c in vis.columns if c != "total"]
    ax.stackplot(vis.index, [vis[c] for c in syss],
                 labels=[SYSTEMS.get(c, c) for c in syss],
                 colors=[colors.get(c) for c in syss], alpha=.85)
    ax.plot(vis.index, vis["total"], "k-", lw=.8, label="Total")
    ax.set(xlabel="Observation time (file scale)", ylabel="Satellites tracked",
           title=f"Satellite visibility over time {title}".strip(),
           ylim=(0, vis["total"].max() * 1.15))
    ax.xaxis.set_major_formatter(fmt)
    ax.grid(alpha=.3)
    ax.legend(loc="upper center", bbox_to_anchor=(.5, -0.12), ncol=4,
              fontsize=8, frameon=False)
    p = outdir / "satellite_count.png"
    fig.tight_layout()
    fig.savefig(p, dpi=150)
    plt.close(fig)
    files.append(p)

    # 4. per-satellite visibility timeline
    v = df.drop_duplicates(["time", "sat"])
    sats = sorted(v["sat"].unique())
    idx = {s_: k for k, s_ in enumerate(sats)}
    fig, ax = plt.subplots(figsize=(11, max(4, .22 * len(sats) + 1.5)))
    ax.scatter(v["time"], v["sat"].map(idx), s=6, marker="|",
               c=[colors.get(x) for x in v["system"]])
    ax.set_yticks(range(len(sats)))
    ax.set_yticklabels(sats, fontsize=7)
    ax.set(xlabel="Observation time (file scale)", title="Satellite tracking timeline")
    ax.xaxis.set_major_formatter(fmt)
    ax.grid(alpha=.3, axis="x")
    p = outdir / "satellite_timeline.png"
    fig.tight_layout()
    fig.savefig(p, dpi=150)
    plt.close(fig)
    files.append(p)
    return files


def export_all(obs_path, outdir, systems="GEC") -> dict:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    header, df = read_obs(obs_path, systems)
    summary = summarize(header, df)
    cp = carrier_phase_and_pseudorange(df)
    df.to_csv(outdir / "observations_long.csv", index=False)
    cp.to_csv(outdir / "phase_pseudorange.csv", index=False)
    observations_wide(df).to_csv(outdir / "observations_wide.csv", index=False)
    snr = snr_table(df)
    snr.to_csv(outdir / "snr_stats.csv", index=False)
    visibility(df).to_csv(outdir / "visibility.csv")
    comp = (df.groupby(["system", "code"]).size().rename("n_obs").reset_index())
    comp.to_csv(outdir / "observation_counts.csv", index=False)
    with open(outdir / "summary.json", "w") as fh:
        json.dump(summary, fh, indent=2, default=str)
    summary["plots"] = [str(p) for p in make_plots(df, outdir, f"({header.marker_name})")]
    return summary


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    o = sub.add_parser("obs", help="parse observation file, export CSV/JSON/plots")
    o.add_argument("file"); o.add_argument("-o", "--out", default="rinex_out")
    o.add_argument("--systems", default="GEC")
    n = sub.add_parser("nav", help="parse navigation file, export ephemeris CSV")
    n.add_argument("file"); n.add_argument("-o", "--out", default="rinex_out")
    n.add_argument("--systems", default="GEC")
    a = ap.parse_args(argv)
    if a.cmd == "obs":
        s = export_all(a.file, a.out, a.systems)
        print(json.dumps({k: s[k] for k in s if k not in ("obs_types",)},
                         indent=2, default=str))
    else:
        nav = read_nav(a.file, a.systems)
        out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
        nav.ephemeris.drop(columns=["raw"], errors="ignore").to_csv(out / "ephemeris.csv", index=False)
        nav.sto.to_csv(out / "sto.csv", index=False)
        nav.ion.to_csv(out / "ion.csv", index=False)
        print(f"RINEX {nav.version}: {len(nav.ephemeris)} EPH, {len(nav.sto)} STO, "
              f"{len(nav.ion)} ION records")
        if not nav.ephemeris.empty:
            print(nav.ephemeris.groupby("system").size().to_string())


if __name__ == "__main__":
    sys.exit(main())
