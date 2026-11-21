"""Reusable dual-frequency STEC/VTEC calculations.

The calculation accepts normalized observation and elevation tables, so parser-
specific file formats can be adapted without coupling TEC processing to a parser.
"""
from __future__ import annotations

from datetime import datetime, timedelta
import math
from pathlib import Path
from typing import Dict, Tuple

import pandas as pd

_TEC_SCALE = 1e16
_EARTH_MEAN_RADIUS = 6371000.0
_GPS_UTC_OFFSET_HISTORY = (
    (datetime(1981, 7, 1), 1),
    (datetime(1982, 7, 1), 2),
    (datetime(1983, 7, 1), 3),
    (datetime(1985, 7, 1), 4),
    (datetime(1988, 1, 1), 5),
    (datetime(1990, 1, 1), 6),
    (datetime(1991, 1, 1), 7),
    (datetime(1992, 7, 1), 8),
    (datetime(1993, 7, 1), 9),
    (datetime(1994, 7, 1), 10),
    (datetime(1996, 1, 1), 11),
    (datetime(1997, 7, 1), 12),
    (datetime(1999, 1, 1), 13),
    (datetime(2006, 1, 1), 14),
    (datetime(2009, 1, 1), 15),
    (datetime(2012, 7, 1), 16),
    (datetime(2015, 7, 1), 17),
    (datetime(2017, 1, 1), 18),
)

OBSERVATION_COLUMNS = (
    "epoch",
    "satellite",
    "code_1_m",
    "code_2_m",
    "phase_1_cycles",
    "phase_2_cycles",
    "frequency_1_hz",
    "frequency_2_hz",
    "lli_1",
    "lli_2",
    "code_pair",
)
ELEVATION_COLUMNS = ("epoch", "satellite", "elevation_deg")
TEC_COLUMNS = (
    "epoch",
    "satellite",
    "code_1_m",
    "code_2_m",
    "code_pair",
    "elevation_deg",
    "stec_tecu_uncalibrated",
    "phase_stec_tecu_uncalibrated",
    "phase_smoothed_stec_tecu_uncalibrated",
    "mapping_factor",
    "vtec_tecu_uncalibrated",
    "phase_smoothed_vtec_tecu_uncalibrated",
)


def normalize_long_observations(
    observations: pd.DataFrame,
    signal_pairs: Dict[str, Tuple[str, str, str, str, float, float]],
) -> pd.DataFrame:
    """Adapt RINEX 3/4 long-form rows to the reusable TEC input schema.

    ``signal_pairs`` maps a constellation letter to
    ``(code_1, code_2, phase_1, phase_2, frequency_1_hz, frequency_2_hz)``.
    Explicit signal selection avoids silently mixing tracking attributes.
    """
    required = {"time", "sat", "code", "value", "lli"}
    missing = required - set(observations.columns)
    if missing:
        raise ValueError(
            "long observations are missing required columns: {}".format(
                ", ".join(sorted(missing))
            )
        )
    normalized = []
    for (epoch, satellite), rows in observations.groupby(
        ["time", "sat"], sort=False
    ):
        system = str(satellite)[0]
        pair = signal_pairs.get(system)
        if pair is None:
            continue
        if len(pair) != 6:
            raise ValueError(
                "signal pair for {} must define four codes and two frequencies".format(
                    system
                )
            )
        code_1, code_2, phase_1, phase_2, frequency_1, frequency_2 = pair
        measurements = {
            row.code: (row.value, row.lli)
            for row in rows.itertuples(index=False)
        }
        normalized.append({
            "epoch": epoch,
            "satellite": satellite,
            "code_1_m": measurements.get(code_1, (None, None))[0],
            "code_2_m": measurements.get(code_2, (None, None))[0],
            "phase_1_cycles": measurements.get(phase_1, (None, None))[0],
            "phase_2_cycles": measurements.get(phase_2, (None, None))[0],
            "frequency_1_hz": frequency_1,
            "frequency_2_hz": frequency_2,
            "lli_1": measurements.get(phase_1, (None, None))[1],
            "lli_2": measurements.get(phase_2, (None, None))[1],
            "code_pair": "{}/{}".format(code_1, code_2),
        })
    return pd.DataFrame.from_records(normalized, columns=OBSERVATION_COLUMNS)


def stec_from_pseudoranges(
    code_1_m: float,
    code_2_m: float,
    frequency_1_hz: float,
    frequency_2_hz: float,
) -> float:
    """Return geometry-free code STEC in TECU."""
    if frequency_1_hz == frequency_2_hz:
        raise ValueError("the two signal frequencies must be different")
    frequency_factor = 40.3 * (
        1 / frequency_2_hz ** 2 - 1 / frequency_1_hz ** 2
    )
    return (code_2_m - code_1_m) / frequency_factor / _TEC_SCALE


def stec_from_carrier_phase(
    phase_1_cycles: float,
    phase_2_cycles: float,
    frequency_1_hz: float,
    frequency_2_hz: float,
) -> float:
    """Return geometry-free carrier-phase STEC in TECU."""
    if frequency_1_hz == frequency_2_hz:
        raise ValueError("the two signal frequencies must be different")
    speed_of_light = 299792458.0
    frequency_factor = 40.3 * (
        1 / frequency_2_hz ** 2 - 1 / frequency_1_hz ** 2
    )
    phase_geometry_free = (
        speed_of_light / frequency_1_hz * phase_1_cycles
        - speed_of_light / frequency_2_hz * phase_2_cycles
    )
    return phase_geometry_free / frequency_factor / _TEC_SCALE


def thin_shell_mapping_factor(elevation_degrees: float, shell_height_km: float) -> float:
    elevation_radians = math.radians(elevation_degrees)
    shell_radius = _EARTH_MEAN_RADIUS + shell_height_km * 1000
    return math.sqrt(
        1 - (_EARTH_MEAN_RADIUS / shell_radius * math.cos(elevation_radians)) ** 2
    )


def observation_time_to_utc(timestamp: datetime, time_system: str) -> datetime:
    """Convert supported GNSS time scales to UTC."""
    scale = time_system.upper()
    if scale in ("UTC", "GLO"):
        return timestamp
    if scale in ("GPS", "GPST", "GAL", "GST"):
        offset_adjustment = 0
    elif scale in ("BDT", "BDS"):
        offset_adjustment = -14
    elif scale == "TAI":
        offset_adjustment = 19
    else:
        raise ValueError(
            "cannot convert observation time system {!r} to UTC".format(time_system)
        )

    utc = timestamp
    for _ in range(3):
        offset = 0
        for effective_date, leap_seconds in _GPS_UTC_OFFSET_HISTORY:
            if utc >= effective_date:
                offset = leap_seconds + offset_adjustment
            else:
                break
        utc = timestamp - timedelta(seconds=offset)
    return utc


def calculate_tec_dataframe(
    observations: pd.DataFrame,
    elevations: pd.DataFrame,
    shell_height_km: float = 350.0,
    time_system: str = "GPS",
    smoothing_window_seconds: float = 100.0,
    code_dcb_corrections_m: Dict[str, float] = None,
) -> pd.DataFrame:
    """Calculate uncalibrated STEC/VTEC from normalized observation tables.

    ``observations`` must contain :data:`OBSERVATION_COLUMNS`; ``elevations``
    must contain :data:`ELEVATION_COLUMNS`. Epochs in both tables must use the
    same time scale. Frequencies are supplied per observation row, allowing
    callers to normalize different GNSS systems and signal pairs.
    """
    if not math.isfinite(shell_height_km) or shell_height_km <= 0:
        raise ValueError("ionosphere shell height must be a positive number")
    if not math.isfinite(smoothing_window_seconds) or smoothing_window_seconds <= 0:
        raise ValueError("smoothing window must be a positive number of seconds")
    missing_observations = set(OBSERVATION_COLUMNS) - set(observations.columns)
    if missing_observations:
        raise ValueError(
            "normalized observations are missing required columns: {}".format(
                ", ".join(sorted(missing_observations))
            )
        )
    missing_elevations = set(ELEVATION_COLUMNS) - set(elevations.columns)
    if missing_elevations:
        raise ValueError(
            "normalized elevations are missing required columns: {}".format(
                ", ".join(sorted(missing_elevations))
            )
        )

    observations = observations.copy()
    code_dcb_corrections_m = code_dcb_corrections_m or {}
    elevations = elevations.loc[:, ELEVATION_COLUMNS].copy()
    observations["epoch"] = pd.to_datetime(observations["epoch"])
    elevations["epoch"] = pd.to_datetime(elevations["epoch"])
    observations = observations.sort_values(["satellite", "epoch"])

    records = []
    for satellite, satellite_rows in observations.groupby("satellite", sort=False):
        sample_intervals = satellite_rows["epoch"].diff().dt.total_seconds()
        sample_intervals = sample_intervals[sample_intervals > 0]
        sample_interval = (
            float(sample_intervals.median()) if not sample_intervals.empty else 1.0
        )
        maximum_count = max(
            1, int(math.ceil(smoothing_window_seconds / sample_interval))
        )
        maximum_gap = max(1.5, 1.5 * sample_interval)
        previous = None
        for row in satellite_rows.itertuples(index=False):
            code_1 = row.code_1_m
            code_2 = row.code_2_m
            phase_1 = row.phase_1_cycles
            phase_2 = row.phase_2_cycles
            frequency_1 = row.frequency_1_hz
            frequency_2 = row.frequency_2_hz
            if (
                pd.isna(code_1)
                or pd.isna(code_2)
                or pd.isna(phase_1)
                or pd.isna(phase_2)
                or code_1 <= 0
                or code_2 <= 0
                or phase_1 <= 0
                or phase_2 <= 0
                or pd.isna(frequency_1)
                or pd.isna(frequency_2)
                or not all(
                    math.isfinite(float(value))
                    for value in (code_1, code_2, phase_1, phase_2,
                                  frequency_1, frequency_2)
                )
                or frequency_1 <= 0
                or frequency_2 <= 0
                or frequency_1 == frequency_2
            ):
                previous = None
                continue

            code_stec = stec_from_pseudoranges(
                code_1,
                code_2 + code_dcb_corrections_m.get(satellite, 0.0),
                frequency_1,
                frequency_2,
            )
            phase_stec = stec_from_carrier_phase(
                phase_1, phase_2, frequency_1, frequency_2
            )
            lli_1 = 0 if pd.isna(row.lli_1) else int(row.lli_1)
            lli_2 = 0 if pd.isna(row.lli_2) else int(row.lli_2)
            lli_slip = (lli_1 & 1) or (lli_2 & 1)
            if (
                previous is None
                or lli_slip
                or (row.epoch - previous["epoch"]).total_seconds() > maximum_gap
            ):
                smoothed_stec = code_stec
                count = 1
            else:
                count = min(previous["count"] + 1, maximum_count)
                alpha = 1.0 / count
                smoothed_stec = alpha * code_stec + (1 - alpha) * (
                    previous["smoothed_stec"] + phase_stec - previous["phase_stec"]
                )
            previous = {
                "epoch": row.epoch,
                "phase_stec": phase_stec,
                "smoothed_stec": smoothed_stec,
                "count": count,
            }
            records.append({
                "epoch": row.epoch,
                "satellite": satellite,
                "code_1_m": code_1,
                "code_2_m": code_2,
                "code_pair": row.code_pair,
                "phase_stec_tecu_uncalibrated": phase_stec,
                "phase_smoothed_stec_tecu_uncalibrated": smoothed_stec,
                "stec_tecu_uncalibrated": code_stec,
            })

    if not records:
        return pd.DataFrame(columns=TEC_COLUMNS)
    calculated = pd.DataFrame.from_records(records)
    result = calculated.merge(
        elevations,
        on=["epoch", "satellite"],
        how="inner",
        validate="many_to_one",
    )
    result = result[result["elevation_deg"].between(0, 90)].copy()
    result["mapping_factor"] = result["elevation_deg"].map(
        lambda value: thin_shell_mapping_factor(value, shell_height_km)
    )
    result["vtec_tecu_uncalibrated"] = (
        result["stec_tecu_uncalibrated"] * result["mapping_factor"]
    )
    result["phase_smoothed_vtec_tecu_uncalibrated"] = (
        result["phase_smoothed_stec_tecu_uncalibrated"] * result["mapping_factor"]
    )
    result["epoch"] = result["epoch"].map(
        lambda value: observation_time_to_utc(value.to_pydatetime(), time_system)
    )
    return result.sort_values(["epoch", "satellite"]).loc[:, TEC_COLUMNS].reset_index(
        drop=True
    )


def elevations_from_compact_file(path) -> Tuple[float, Dict[datetime, Dict[str, float]]]:
    """Read TEQC-style COMPACT2 elevation rows into epoch-indexed dictionaries."""
    path = Path(path)
    with path.open("r", encoding="ascii") as stream:
        format_line = stream.readline().strip()
        interval_line = stream.readline().split()
        start_line = stream.readline().split()
        if format_line != "COMPACT2":
            raise ValueError("{} is not a COMPACT2 elevation file".format(path))
        if len(interval_line) != 2 or interval_line[0] != "T_SAMP":
            raise ValueError("elevation file has an invalid T_SAMP header")
        if len(start_line) != 2 or start_line[0] != "START_TIME_MJD":
            raise ValueError("elevation file has an invalid START_TIME_MJD header")
        try:
            interval = float(interval_line[1])
            start_mjd = float(start_line[1])
        except ValueError as error:
            raise ValueError("elevation file has invalid time metadata") from error
        if interval <= 0:
            raise ValueError("elevation sampling interval must be positive")

        start_time = datetime(1858, 11, 17) + timedelta(days=start_mjd)
        start_time = (start_time + timedelta(microseconds=500000)).replace(
            microsecond=0
        )
        elevations = {}
        satellites = None
        epoch_index = 0
        while True:
            satellite_line = stream.readline()
            if not satellite_line:
                break
            angle_line = stream.readline()
            if not angle_line:
                raise ValueError("incomplete elevation row at end of file")
            marker = satellite_line.split()
            if not marker:
                raise ValueError("empty satellite list in elevation file")
            if marker[0] != "-1":
                try:
                    count = int(marker[0])
                except ValueError as error:
                    raise ValueError("invalid satellite count in elevation file") from error
                satellites = marker[1:]
                if len(satellites) != count:
                    raise ValueError("satellite count does not match elevation list")
            elif satellites is None:
                raise ValueError("elevation file uses -1 before declaring satellites")
            try:
                values = [float(value) for value in angle_line.split()]
            except ValueError as error:
                raise ValueError("invalid elevation value") from error
            if len(values) != len(satellites):
                raise ValueError("elevation count does not match satellite list")
            timestamp = start_time + timedelta(seconds=epoch_index * interval)
            elevations[timestamp] = dict(zip(satellites, values))
            epoch_index += 1

    return interval, elevations
