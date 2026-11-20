#!/usr/bin/env python3
"""
rinex_tool - RINEX 3.x / 4.x observation and navigation parser, time-scale
conversion, tidy pandas extraction and summary statistics.

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
from datetime import date, timedelta
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
from ..tec_calculation import (
    calculate_tec_dataframe as calculate_normalized_tec_dataframe,
    elevations_from_compact_file,
    normalize_long_observations,
)

SYSTEMS = {"G": "GPS", "E": "Galileo", "C": "BeiDou", "R": "GLONASS",
           "J": "QZSS", "I": "NavIC", "S": "SBAS"}
_GPS_EPOCH = pd.Timestamp("1980-01-06")
_GPS_WEEK_SECONDS = 604800.0
_SPEED_OF_LIGHT = 299792458.0
_EARTH_GRAVITATIONAL_CONSTANT = {
    "G": 3.986005e14,
    "E": 3.986004418e14,
    "C": 3.986004418e14,
}
_EARTH_ROTATION_RATE = 7.2921151467e-5
_SUPPORTED_ORBIT_SYSTEMS = frozenset(_EARTH_GRAVITATIONAL_CONSTANT)


def _version_output_dir(outdir, version: float) -> Path:
    return Path(outdir) / "v{:03d}".format(round(version * 100))


def _project_data_dir() -> Path:
    root = Path(__file__).resolve().parents[3]
    return root if (root / "Data_3_4").exists() else Path.cwd()


def _dated_path(outdir: Path, filename: str, date_stamp: str) -> Path:
    path = Path(filename)
    return outdir / "{}_{}{}".format(path.stem, date_stamp, path.suffix)


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
            tokens = line[:60].split()
            if len(tokens) >= 5:
                code, a0, a1, reference_time, week = tokens[:5]
            else:
                code = line[:4].strip()
                a0, a1, reference_time, week = (
                    line[4:22],
                    line[22:38],
                    line[38:45],
                    line[45:50],
                )
            if code and str(a0).strip() and str(a1).strip() and str(reference_time).strip():
                hdr["time_corr"][code] = {
                    "a0": _f(a0),
                    "a1": _f(a1),
                    "t": int(float(reference_time)),
                    "w": int(float(week)) if str(week).strip() else 0,
                }
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


# ------------------------------------------------------------------ statistics
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


def _satellite_position_gnss(
    timestamp_gpst: pd.Timestamp,
    satellite: str,
    ephemerides: pd.DataFrame,
    pseudorange: Optional[float],
    ephemeris_records=None,
) -> tuple:
    system = satellite[0]
    if system not in _SUPPORTED_ORBIT_SYSTEMS:
        raise ValueError("broadcast orbit calculation is not supported for {}".format(system))
    candidates = []
    if ephemeris_records is None:
        records = ephemerides[ephemerides["sat"] == satellite]
        if records.empty:
            raise ValueError("navigation file has no ephemeris for {}".format(satellite))
        record_rows = []
        required = ("sqrt_a", "e", "m0", "delta_n", "omega", "omega0", "i0",
                    "toe", "toc_gpst")
        for _, eph in records.iterrows():
            if pd.isna(eph.get("sv_health", 0)) or float(eph.get("sv_health", 0)) != 0:
                continue
            if any(pd.isna(eph.get(name)) for name in required):
                continue
            toc_gpst = pd.Timestamp(eph["toc_gpst"])
            week_start = _GPS_EPOCH + pd.to_timedelta(
                np.floor((toc_gpst - _GPS_EPOCH).total_seconds() / _GPS_WEEK_SECONDS)
                * _GPS_WEEK_SECONDS,
                unit="s",
            )
            toe_seconds = float(eph["toe"]) % _GPS_WEEK_SECONDS
            system_time_offset = 14.0 if system == "C" else 0.0
            toe_gpst = min(
                (
                    week_start + pd.to_timedelta(
                        toe_seconds + system_time_offset + offset * _GPS_WEEK_SECONDS,
                        unit="s",
                    )
                    for offset in (-1, 0, 1)
                ),
                key=lambda value: abs(value - toc_gpst),
            )
            record_rows.append((toe_gpst, eph))
    else:
        record_rows = ephemeris_records

    for toe_gpst, eph in record_rows:
        tk = (timestamp_gpst - toe_gpst).total_seconds()
        tk = (tk + _GPS_WEEK_SECONDS / 2) % _GPS_WEEK_SECONDS - _GPS_WEEK_SECONDS / 2
        candidates.append((abs(tk), tk, eph))
    if not candidates:
        raise ValueError("navigation file has no healthy orbital ephemeris for {}".format(satellite))

    _, tk, eph = min(candidates, key=lambda value: value[0])
    sqrt_a = float(eph["sqrt_a"])
    eccentricity = float(eph["e"])
    if sqrt_a <= 0 or not 0 <= eccentricity < 1:
        raise ValueError("invalid orbital parameters for {}".format(satellite))
    semi_major_axis = sqrt_a * sqrt_a
    mean_motion = np.sqrt(_EARTH_GRAVITATIONAL_CONSTANT[system] / semi_major_axis ** 3)
    mean_anomaly = float(eph["m0"]) + (
        mean_motion + float(eph.get("delta_n", 0.0))
    ) * tk
    eccentric_anomaly = mean_anomaly
    for _ in range(20):
        updated = mean_anomaly + eccentricity * np.sin(eccentric_anomaly)
        if abs(updated - eccentric_anomaly) < 1e-13:
            eccentric_anomaly = updated
            break
        eccentric_anomaly = updated

    true_anomaly = np.arctan2(
        np.sqrt(1 - eccentricity * eccentricity) * np.sin(eccentric_anomaly),
        np.cos(eccentric_anomaly) - eccentricity,
    )
    argument = true_anomaly + float(eph["omega"])
    twice_argument = 2 * argument
    argument += (
        float(eph.get("cus", 0.0)) * np.sin(twice_argument)
        + float(eph.get("cuc", 0.0)) * np.cos(twice_argument)
    )
    radius = (
        semi_major_axis * (1 - eccentricity * np.cos(eccentric_anomaly))
        + float(eph.get("crs", 0.0)) * np.sin(twice_argument)
        + float(eph.get("crc", 0.0)) * np.cos(twice_argument)
    )
    inclination = (
        float(eph["i0"])
        + float(eph.get("idot", 0.0)) * tk
        + float(eph.get("cis", 0.0)) * np.sin(twice_argument)
        + float(eph.get("cic", 0.0)) * np.cos(twice_argument)
    )
    orbital_x = radius * np.cos(argument)
    orbital_y = radius * np.sin(argument)
    longitude = (
        float(eph["omega0"])
        + (float(eph.get("omega_dot", 0.0)) - _EARTH_ROTATION_RATE) * tk
        - _EARTH_ROTATION_RATE * float(eph["toe"])
    )
    x = orbital_x * np.cos(longitude) - orbital_y * np.cos(inclination) * np.sin(longitude)
    y = orbital_x * np.sin(longitude) + orbital_y * np.cos(inclination) * np.cos(longitude)
    z = orbital_y * np.sin(inclination)

    if system == "C" and (1 <= int(satellite[1:]) <= 5 or int(satellite[1:]) >= 59):
        rotation = _EARTH_ROTATION_RATE * tk
        tilt = np.radians(-5.0)
        x, y, z = (
            x * np.cos(rotation) + y * np.sin(rotation) * np.cos(tilt)
            + z * np.sin(rotation) * np.sin(tilt),
            -x * np.sin(rotation) + y * np.cos(rotation) * np.cos(tilt)
            + z * np.cos(rotation) * np.sin(tilt),
            -y * np.sin(tilt) + z * np.cos(tilt),
        )

    if pseudorange is not None and pseudorange > 0:
        rotation = _EARTH_ROTATION_RATE * pseudorange / _SPEED_OF_LIGHT
        x, y = (
            x * np.cos(rotation) + y * np.sin(rotation),
            -x * np.sin(rotation) + y * np.cos(rotation),
        )
    return float(x), float(y), float(z)


def _azimuth_elevation_gnss(satellite_xyz: tuple, receiver_xyz: tuple) -> tuple:
    x, y, z = receiver_xyz
    longitude = np.arctan2(y, x)
    latitude = np.arctan2(z, np.hypot(x, y) * (1 - 6.69437999014e-3))
    dx, dy, dz = (satellite_xyz[index] - receiver_xyz[index] for index in range(3))
    east = -np.sin(longitude) * dx + np.cos(longitude) * dy
    north = (
        -np.sin(latitude) * np.cos(longitude) * dx
        - np.sin(latitude) * np.sin(longitude) * dy
        + np.cos(latitude) * dz
    )
    up = (
        np.cos(latitude) * np.cos(longitude) * dx
        + np.cos(latitude) * np.sin(longitude) * dy
        + np.sin(latitude) * dz
    )
    azimuth = np.degrees(np.arctan2(east, north))
    elevation = np.degrees(np.arctan2(up, np.hypot(east, north)))
    return float(azimuth), float(elevation)


def write_azimuth_elevation_files(
    observation_path,
    navigation_path,
    output_dir=None,
    sample_interval: int = 15,
    systems: str = "GEC",
) -> tuple[Path, Path]:
    """Write compact azimuth/elevation files from RINEX 3/4 observation/nav data."""
    if sample_interval <= 0:
        raise ValueError("sample interval must be a positive number of seconds")
    observation_path = Path(observation_path)
    navigation_path = Path(navigation_path)
    output_root = (
        Path(output_dir)
        if output_dir is not None
        else _project_data_dir() / "Data_3_4" / "output"
    )
    header, observations = read_obs(observation_path, systems)
    if header.approx_xyz is None:
        raise ValueError("observation header is missing APPROX POSITION XYZ")
    nav = read_nav(navigation_path, systems)
    ephemerides = nav.ephemeris
    if ephemerides.empty:
        raise ValueError("navigation file contains no supported ephemerides")

    output_dir = _version_output_dir(output_root, header.version)
    output_dir.mkdir(parents=True, exist_ok=True)
    azimuth_path = output_dir / (observation_path.stem + ".azi")
    elevation_path = output_dir / (observation_path.stem + ".ele")
    epochs = pd.DatetimeIndex(observations["time"].drop_duplicates().sort_values())
    if epochs.empty:
        raise ValueError("observation file contains no usable epochs")
    first_epoch = epochs[0]
    previous_satellites = None
    receiver_xyz = header.approx_xyz
    navigation_satellites = _supported_satellites(ephemerides, systems)
    if not navigation_satellites:
        raise ValueError("navigation file has no healthy GPS/Galileo/BeiDou ephemerides")
    ephemerides_by_satellite = _prepare_ephemerides(ephemerides)
    satellites_by_epoch = (
        observations.loc[:, ["time", "sat"]]
        .drop_duplicates()
        .groupby("time", sort=False)["sat"]
        .agg(tuple)
        .to_dict()
    )
    pseudorange_rows = observations.loc[
        (observations["type"] == "C") & (observations["value"] > 0),
        ["time", "sat", "value"],
    ].drop_duplicates(["time", "sat"])
    pseudoranges_by_epoch = {
        epoch: dict(zip(rows["sat"], rows["value"]))
        for epoch, rows in pseudorange_rows.groupby("time", sort=False)
    }

    with (
        azimuth_path.open("w", encoding="ascii", newline="\n") as azimuth_file,
        elevation_path.open("w", encoding="ascii", newline="\n") as elevation_file,
    ):
        for output_file in (azimuth_file, elevation_file):
            output_file.write("COMPACT2\n")
            output_file.write("T_SAMP    {:0.4f}\n".format(sample_interval))
            mjd = (
                first_epoch - pd.Timestamp("1858-11-17")
            ).total_seconds() / 86400.0
            output_file.write("START_TIME_MJD  {:0.9f}\n".format(mjd))

        for epoch in epochs:
            elapsed = (epoch - first_epoch).total_seconds()
            if round(elapsed) % sample_interval:
                continue
            observed_satellites = set(satellites_by_epoch.get(epoch, ()))
            satellites_to_process = sorted(
                observed_satellites.intersection(navigation_satellites)
            )
            pseudoranges = pseudoranges_by_epoch.get(epoch, {})
            epoch_gpst = to_gpst(pd.DatetimeIndex([epoch]), header.time_system)[0]
            azimuths = []
            elevations = []
            for satellite in satellites_to_process:
                position = _satellite_position_gnss(
                    epoch_gpst,
                    satellite,
                    ephemerides,
                    pseudoranges.get(satellite),
                    ephemerides_by_satellite.get(satellite),
                )
                azimuth, elevation = _azimuth_elevation_gnss(position, receiver_xyz)
                if elevation < 0:
                    continue
                azimuths.append((satellite, azimuth))
                elevations.append((satellite, elevation))

            satellites = tuple(satellite for satellite, _ in azimuths)
            marker = (
                "{} {}".format(len(satellites), " ".join(satellites))
                if satellites != previous_satellites
                else "-1"
            )
            azimuth_file.write(marker + "\n")
            elevation_file.write(marker + "\n")
            azimuth_file.write(
                "".join("{:8.3f} ".format(value) for _, value in azimuths) + "\n"
            )
            elevation_file.write(
                "".join("{:8.3f} ".format(value) for _, value in elevations) + "\n"
            )
            previous_satellites = satellites
    return azimuth_path, elevation_path


def calculate_tec_from_observation_file(
    observation_path,
    elevation_path,
    system: str = "G",
    code_1: str = "C1C",
    code_2: str = "C2W",
    phase_1: str = "L1C",
    phase_2: str = "L2W",
    frequency_1_hz: float = 1575.42e6,
    frequency_2_hz: float = 1227.60e6,
    shell_height_km: float = 350.0,
    smoothing_window_seconds: float = 100.0,
) -> pd.DataFrame:
    """Calculate TEC from one explicitly selected constellation signal pair."""
    header, observations = read_obs(observation_path, systems=system)
    if system not in header.obs_types:
        raise ValueError(
            "observation file does not contain constellation {}".format(system)
        )
    missing_signals = {
        code_1, code_2, phase_1, phase_2
    } - set(header.obs_types[system])
    if missing_signals:
        raise ValueError(
            "constellation {} is missing selected signal(s): {}".format(
                system, ", ".join(sorted(missing_signals))
            )
        )
    _, elevations_by_epoch = elevations_from_compact_file(elevation_path)
    elevation_rows = [
        {"epoch": timestamp, "satellite": satellite, "elevation_deg": elevation}
        for timestamp, satellites in elevations_by_epoch.items()
        for satellite, elevation in satellites.items()
    ]
    normalized = normalize_long_observations(
        observations,
        {
            system: (
                code_1,
                code_2,
                phase_1,
                phase_2,
                frequency_1_hz,
                frequency_2_hz,
            )
        },
    )
    return calculate_normalized_tec_dataframe(
        normalized,
        pd.DataFrame.from_records(
            elevation_rows, columns=("epoch", "satellite", "elevation_deg")
        ),
        shell_height_km=shell_height_km,
        time_system=header.time_system,
        smoothing_window_seconds=smoothing_window_seconds,
    )


def _supported_satellites(ephemerides: pd.DataFrame, systems: str):
    allowed = set(systems) & _SUPPORTED_ORBIT_SYSTEMS
    if ephemerides.empty:
        return set()
    available = ephemerides[ephemerides["system"].isin(allowed)]
    if "sv_health" in available:
        available = available[available["sv_health"].fillna(1) == 0]
    return set(available["sat"])


def _prepare_ephemerides(ephemerides: pd.DataFrame) -> dict:
    required = ("sqrt_a", "e", "m0", "delta_n", "omega", "omega0", "i0",
                "toe", "toc_gpst")
    prepared = {}
    for satellite, records in ephemerides.groupby("sat", sort=False):
        system_time_offset = 14.0 if satellite[0] == "C" else 0.0
        prepared_records = []
        for eph in records.to_dict("records"):
            if (
                pd.isna(eph.get("sv_health", 0))
                or float(eph.get("sv_health", 0)) != 0
                or any(pd.isna(eph.get(name)) for name in required)
            ):
                continue
            toc_gpst = pd.Timestamp(eph["toc_gpst"])
            week_number = np.floor(
                (toc_gpst - _GPS_EPOCH).total_seconds() / _GPS_WEEK_SECONDS
            )
            week_start = _GPS_EPOCH + pd.to_timedelta(
                week_number * _GPS_WEEK_SECONDS, unit="s"
            )
            toe_seconds = float(eph["toe"]) % _GPS_WEEK_SECONDS
            toe_gpst = min(
                (
                    week_start + pd.to_timedelta(
                        toe_seconds + system_time_offset + offset * _GPS_WEEK_SECONDS,
                        unit="s",
                    )
                    for offset in (-1, 0, 1)
                ),
                key=lambda value: abs(value - toc_gpst),
            )
            prepared_records.append((toe_gpst, eph))
        if prepared_records:
            prepared[satellite] = prepared_records
    return prepared


def export_all(obs_path, outdir, systems="GEC") -> dict:
    header, df = read_obs(obs_path, systems)
    outdir = _version_output_dir(outdir, header.version)
    outdir.mkdir(parents=True, exist_ok=True)
    date_stamp = date.today().strftime("%d%m%y")
    summary = summarize(header, df)
    cp = carrier_phase_and_pseudorange(df)
    df.to_csv(_dated_path(outdir, "observations_long.csv", date_stamp), index=False)
    cp.to_csv(_dated_path(outdir, "phase_pseudorange.csv", date_stamp), index=False)
    observations_wide(df).to_csv(
        _dated_path(outdir, "observations_wide.csv", date_stamp), index=False
    )
    snr = snr_table(df)
    snr.to_csv(_dated_path(outdir, "snr_stats.csv", date_stamp), index=False)
    visibility(df).to_csv(_dated_path(outdir, "visibility.csv", date_stamp))
    comp = (df.groupby(["system", "code"]).size().rename("n_obs").reset_index())
    comp.to_csv(_dated_path(outdir, "observation_counts.csv", date_stamp), index=False)
    with open(_dated_path(outdir, "summary.json", date_stamp), "w") as fh:
        json.dump(summary, fh, indent=2, default=str)
    return summary


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    o = sub.add_parser("obs", help="parse observation file and export CSV/JSON")
    default_out = _project_data_dir() / "Data_3_4" / "output"
    o.add_argument("file"); o.add_argument("-o", "--out", default=default_out)
    o.add_argument("--systems", default="GEC")
    n = sub.add_parser("nav", help="parse navigation file, export ephemeris CSV")
    n.add_argument("file"); n.add_argument("-o", "--out", default=default_out)
    n.add_argument("--systems", default="GEC")
    g = sub.add_parser(
        "angles", help="generate compact azimuth/elevation from observation and nav files"
    )
    g.add_argument("obs_file")
    g.add_argument("nav_file")
    g.add_argument("-o", "--out", default=default_out)
    g.add_argument("--systems", default="GEC")
    g.add_argument("--interval", type=int, default=15)
    t = sub.add_parser(
        "tec", help="calculate STEC/VTEC from observations and compact elevation"
    )
    t.add_argument("obs_file")
    t.add_argument("--ele", type=Path)
    t.add_argument("--out", type=Path)
    t.add_argument("--system", choices=("G", "E", "C"), default="G")
    t.add_argument("--code1")
    t.add_argument("--code2")
    t.add_argument("--phase1")
    t.add_argument("--phase2")
    t.add_argument("--frequency1", type=float, metavar="HZ")
    t.add_argument("--frequency2", type=float, metavar="HZ")
    t.add_argument("--shell-height", type=float, default=350.0, metavar="KM")
    t.add_argument("--smoothing-window", type=float, default=100.0, metavar="SECONDS")
    a = ap.parse_args(argv)
    if a.cmd == "obs":
        s = export_all(a.file, a.out, a.systems)
        print(json.dumps({k: s[k] for k in s if k not in ("obs_types",)},
                         indent=2, default=str))
    elif a.cmd == "angles":
        try:
            azimuth_path, elevation_path = write_azimuth_elevation_files(
                a.obs_file, a.nav_file, a.out, a.interval, a.systems
            )
        except (OSError, ValueError) as error:
            print("Error: {}".format(error), file=sys.stderr)
            return 2
        print("Azimuth written to: {}".format(azimuth_path))
        print("Elevation written to: {}".format(elevation_path))
    elif a.cmd == "tec":
        observation_path = Path(a.obs_file)
        try:
            header, _ = read_obs(observation_path, systems=a.system)
            signal_defaults = {
                "G": ("C1C", "C2W", "L1C", "L2W", 1575.42e6, 1227.60e6),
                "E": ("C1C", "C5Q", "L1C", "L5Q", 1575.42e6, 1176.45e6),
                "C": ("C2I", "C7I", "L2I", "L7I", 1561.098e6, 1207.14e6),
            }
            defaults = signal_defaults[a.system]
            selected = (
                a.code1 or defaults[0],
                a.code2 or defaults[1],
                a.phase1 or defaults[2],
                a.phase2 or defaults[3],
                a.frequency1 or defaults[4],
                a.frequency2 or defaults[5],
            )
            explicit_csv = a.out and Path(a.out).suffix.lower() == ".csv"
            output_dir = (
                Path(a.out).parent
                if explicit_csv
                else _version_output_dir(a.out or default_out, header.version)
            )
            elevation_path = a.ele or (output_dir / (observation_path.stem + ".ele"))
            output_path = (
                Path(a.out)
                if explicit_csv
                else output_dir / (observation_path.stem + "_tec.csv")
            )
            tec = calculate_tec_from_observation_file(
                observation_path,
                elevation_path,
                system=a.system,
                code_1=selected[0],
                code_2=selected[1],
                phase_1=selected[2],
                phase_2=selected[3],
                frequency_1_hz=selected[4],
                frequency_2_hz=selected[5],
                shell_height_km=a.shell_height,
                smoothing_window_seconds=a.smoothing_window,
            )
            if tec.empty:
                raise ValueError(
                    "no matching epochs with dual-frequency observations and elevation"
                )
            output_path.parent.mkdir(parents=True, exist_ok=True)
            tec.to_csv(output_path, index=False)
        except (OSError, ValueError, KeyError) as error:
            print("Error: {}".format(error), file=sys.stderr)
            return 2
        print("TEC rows: {}".format(len(tec)))
        print("Output written to: {}".format(output_path))
        print(
            "Generate plots with: python3 rinex_visualization.py tec {}".format(
                output_path
            )
        )
    else:
        nav = read_nav(a.file, a.systems)
        out = _version_output_dir(a.out, nav.version)
        out.mkdir(parents=True, exist_ok=True)
        date_stamp = date.today().strftime("%d%m%y")
        nav.ephemeris.drop(columns=["raw"], errors="ignore").to_csv(
            _dated_path(out, "ephemeris.csv", date_stamp), index=False
        )
        nav.sto.to_csv(_dated_path(out, "sto.csv", date_stamp), index=False)
        nav.ion.to_csv(_dated_path(out, "ion.csv", date_stamp), index=False)
        print(f"RINEX {nav.version}: {len(nav.ephemeris)} EPH, {len(nav.sto)} STO, "
              f"{len(nav.ion)} ION records")
        if not nav.ephemeris.empty:
            print(nav.ephemeris.groupby("system").size().to_string())
    return 0


if __name__ == "__main__":
    sys.exit(main())
