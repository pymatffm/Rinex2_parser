#!/usr/bin/env python3
"""Parse and summarize a RINEX 2 observation file."""

import argparse
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime, timedelta
import math
from pathlib import Path
import sys
from typing import Dict, Iterator, List, Optional, TextIO, Tuple

import numpy as np
import pandas as pd
from tec_calculation import (
    ELEVATION_COLUMNS,
    OBSERVATION_COLUMNS,
    calculate_tec_dataframe as calculate_normalized_tec_dataframe,
    elevations_from_compact_file,
)


NAVIGATION_FIELDS = (
    "af0",
    "af1",
    "af2",
    "iode",
    "crs",
    "delta_n",
    "m0",
    "cuc",
    "ecc",
    "cus",
    "sqrt_a",
    "toe",
    "cic",
    "omega0",
    "cis",
    "inc",
    "crc",
    "omega",
    "dot_omega",
    "idot",
    "codes",
    "gps_week",
    "l2p",
    "sv_accuracy",
    "sv_health",
    "tgd",
    "iodc",
    "transmit_time",
    "fit_interval",
)
_GPS_EPOCH = datetime(1980, 1, 6)
_GPS_WEEK_SECONDS = 604800.0
_SPEED_OF_LIGHT = 299792458.0
_EARTH_GRAVITATIONAL_CONSTANT = 3.986005e14
_EARTH_ROTATION_RATE = 7.2921151467e-5
_RELATIVISTIC_CORRECTION = -4.442807633e-10

_HEADER_LABELS = {
    "ANT # / TYPE",
    "ANTENNA: DELTA H/E/N",
    "APPROX POSITION XYZ",
    "COMMENT",
    "LEAP SECONDS",
    "MARKER NAME",
    "MARKER NUMBER",
    "OBSERVER / AGENCY",
    "REC # / TYPE / VERS",
    "TIME OF FIRST OBS",
    "TIME OF LAST OBS",
    "WAVELENGTH FACT L1/2",
    "# / TYPES OF OBSERV",
}

DATA_2_DIR = Path(__file__).resolve().parent / "Data_2"
DATA_2_OUTPUT_DIR = DATA_2_DIR / "output"


@dataclass(frozen=True)
class RinexHeader:
    version: float
    file_type: str
    marker_name: str
    observation_types: Tuple[str, ...]
    approximate_position_xyz: Optional[Tuple[float, float, float]] = None
    time_system: str = "GPS"


@dataclass(frozen=True)
class Observation:
    value: Optional[float]
    loss_of_lock_indicator: Optional[int]
    signal_strength: Optional[int]


@dataclass(frozen=True)
class Epoch:
    timestamp: datetime
    flag: int
    observations: Dict[str, Tuple[Observation, ...]]


@dataclass(frozen=True)
class NavigationHeader:
    version: float
    file_type: str
    leap_seconds: Optional[int] = None


@dataclass(frozen=True)
class Ephemeris:
    satellite: str
    toc: datetime
    parameters: Tuple[float, ...]


def _read_header(stream) -> RinexHeader:
    version_line = stream.readline()
    if not version_line:
        raise ValueError("the file is empty")

    try:
        version = float(version_line[:9])
    except ValueError as error:
        raise ValueError("the first line does not contain a RINEX version") from error
    if not 2.0 <= version < 3.0:
        raise ValueError("only RINEX 2 observation files are supported")

    file_type = version_line[20:40].strip()
    marker_name = ""
    approximate_position_xyz = None
    time_system = "GPS"
    observation_types = []
    observation_type_count = None

    while True:
        line = stream.readline()
        if not line:
            raise ValueError("missing END OF HEADER")
        label = line[60:80].strip()

        if label == "MARKER NAME":
            marker_name = line[:60].strip()
        elif label == "APPROX POSITION XYZ":
            try:
                coordinates = tuple(float(value) for value in line[:60].split()[:3])
                if len(coordinates) != 3:
                    raise ValueError
                approximate_position_xyz = coordinates
            except ValueError as error:
                raise ValueError("invalid approximate receiver position") from error
        elif label == "TIME OF FIRST OBS":
            time_system = line[48:51].strip() or "GPS"
        elif label == "# / TYPES OF OBSERV":
            if observation_type_count is None:
                try:
                    observation_type_count = int(line[:6])
                except ValueError as error:
                    raise ValueError("invalid observation-type count in header") from error
            observation_types.extend(
                field.strip() for field in
                (line[start:start + 6] for start in range(6, 60, 6))
                if field.strip()
            )
        elif label == "END OF HEADER":
            break

    if observation_type_count is None:
        raise ValueError("missing # / TYPES OF OBSERV header record")
    observation_types = observation_types[:observation_type_count]
    if len(observation_types) != observation_type_count:
        raise ValueError("header contains fewer observation types than declared")

    return RinexHeader(
        version=version,
        file_type=file_type,
        marker_name=marker_name,
        observation_types=tuple(observation_types),
        approximate_position_xyz=approximate_position_xyz,
        time_system=time_system,
    )


def _read_navigation_header(stream: TextIO) -> NavigationHeader:
    version_line = stream.readline()
    if not version_line:
        raise ValueError("the file is empty")

    try:
        version = float(version_line[:9])
    except ValueError as error:
        raise ValueError("the first line does not contain a RINEX version") from error
    if not 2.0 <= version < 3.0:
        raise ValueError("only RINEX 2 navigation files are supported")

    file_type = version_line[20:40].strip()
    if file_type and "NAVIGATION" not in file_type:
        raise ValueError("not a RINEX 2 navigation file")
    leap_seconds = None
    while True:
        line = stream.readline()
        if not line:
            raise ValueError("missing END OF HEADER")
        label = line[60:80].strip()
        if label == "LEAP SECONDS":
            try:
                leap_seconds = int(line[:6])
            except ValueError as error:
                raise ValueError("invalid leap-second count in header") from error
        elif label == "END OF HEADER":
            break

    return NavigationHeader(
        version=version,
        file_type=file_type,
        leap_seconds=leap_seconds,
    )


def _satellite_id(field: str) -> Optional[str]:
    token = field.strip().upper()
    if not token:
        return None
    if token[0].isalpha():
        return token
    return "G{:02d}".format(int(token))


def _read_satellites(first_line: str, count: int, stream) -> Tuple[str, ...]:
    satellites = []
    line = first_line

    while len(satellites) < count:
        for start in range(32, 68, 3):
            satellite = _satellite_id(line[start:start + 3])
            if satellite is not None:
                satellites.append(satellite)
                if len(satellites) == count:
                    break
        if len(satellites) < count:
            line = stream.readline()
            if not line:
                raise ValueError("unexpected end of file in satellite list")

    return tuple(satellites)


def _read_observation(line: str) -> Observation:
    field = line.ljust(16)
    value_text = field[:14].strip()
    lli_text = field[14].strip()
    strength_text = field[15].strip()

    try:
        value = float(value_text) if value_text else None
        lli = int(lli_text) if lli_text else None
        strength = int(strength_text) if strength_text else None
    except ValueError as error:
        raise ValueError("invalid numeric observation field {!r}".format(field)) from error

    if value is not None and not math.isfinite(value):
        raise ValueError("observation values must be finite")
    return Observation(value, lli, strength)


def _skip_splice_metadata(stream) -> Optional[str]:
    position = stream.tell()
    marker = stream.readline()
    if marker[:60].strip() != "RINEX FILE SPLICE" or marker[60:80].strip() != "COMMENT":
        stream.seek(position)
        return None

    while True:
        line = stream.readline()
        if not line or line[60:80].strip() != "COMMENT":
            return line


def _parse_float_field(field: str) -> float:
    text = field.strip().replace("D", "E").replace("d", "e")
    if not text:
        raise ValueError("missing navigation parameter")
    try:
        value = float(text)
    except ValueError as error:
        raise ValueError("invalid navigation parameter {!r}".format(text)) from error
    if not math.isfinite(value):
        raise ValueError("navigation parameters must be finite")
    return value


def _iter_navigation_records(stream: TextIO) -> Iterator[List[str]]:
    record = []
    for line in stream:
        if line.strip():
            tokens = line[:23].split()
            is_epoch_line = False
            if len(tokens) >= 6:
                try:
                    year = int(tokens[1])
                    month = int(tokens[2])
                    day = int(tokens[3])
                    hour = int(tokens[4])
                    minute = int(tokens[5])
                    is_epoch_line = (
                        0 <= year <= 99
                        and 1 <= month <= 12
                        and 1 <= day <= 31
                        and 0 <= hour <= 23
                        and 0 <= minute <= 59
                    )
                except (ValueError, IndexError):
                    pass
            if is_epoch_line and record:
                yield record
                record = []
            record.append(line.rstrip("\r\n"))
    if record:
        yield record


def _parse_navigation_record(record: List[str]) -> Ephemeris:
    if len(record) != 8:
        raise ValueError(
            "invalid RINEX 2 navigation record: expected 8 lines, got {}".format(
                len(record)
            )
        )

    first = record[0]
    satellite = _satellite_id(first[:3])
    if satellite is None:
        raise ValueError("invalid navigation satellite identifier")

    try:
        if first[0].isalpha():
            clock_start = 23
            date_fields = first[3:23].split()
            year, month, day, hour, minute = map(int, date_fields[:5])
            seconds = float(date_fields[5])
        else:
            clock_start = 22
            year = int(first[2:5])
            month = int(first[5:8])
            day = int(first[8:11])
            hour = int(first[11:14])
            minute = int(first[14:17])
            seconds = float(first[17:22])
        if year < 100:
            year += 1900 if year >= 80 else 2000
        toc = datetime(year, month, day, hour, minute) + timedelta(seconds=seconds)
    except (ValueError, OverflowError) as error:
        raise ValueError("invalid navigation time record") from error

    fields = [
        first[clock_start + start:clock_start + start + 19]
        for start in (0, 19, 38)
    ]
    for line_number, line in enumerate(record[1:], start=1):
        width = 2 if line_number == 7 else 4
        fields.extend(
            line[3 + start:3 + start + 19]
            for start in range(0, width * 19, 19)
        )
    fields = fields[:len(NAVIGATION_FIELDS)]
    if len(fields) != len(NAVIGATION_FIELDS):
        raise ValueError("navigation record contains missing parameters")

    parameters = tuple(_parse_float_field(field) for field in fields)
    return Ephemeris(satellite=satellite, toc=toc, parameters=parameters)


@contextmanager
def open_navigation_file(path: Path):
    """Open a RINEX 2 navigation file and stream ephemeris records."""
    with path.open("r", encoding="ascii", errors="replace", newline="") as stream:
        header = _read_navigation_header(stream)
        yield header, (
            _parse_navigation_record(record)
            for record in _iter_navigation_records(stream)
        )


def read_navigation_dataframe(path: Path):
    """Return the RINEX 2 navigation header and ephemeris DataFrame."""
    records = []
    with open_navigation_file(path) as (header, ephemerides):
        for ephemeris in ephemerides:
            record = dict(zip(NAVIGATION_FIELDS, ephemeris.parameters))
            record["satellite"] = ephemeris.satellite
            record["toc"] = ephemeris.toc
            records.append(record)

    if not records:
        dataframe = pd.DataFrame(
            columns=("satellite", "toc", *NAVIGATION_FIELDS)
        )
    else:
        dataframe = pd.DataFrame(records)
        dataframe = dataframe.set_index(["satellite", "toc"], drop=True)
    return header, dataframe


def _wrap_gps_week(seconds: float) -> float:
    return (seconds + _GPS_WEEK_SECONDS / 2) % _GPS_WEEK_SECONDS - _GPS_WEEK_SECONDS / 2


def _satellite_position(
    timestamp: datetime,
    satellite: str,
    navigation: pd.DataFrame,
    pseudorange: Optional[float],
) -> Tuple[float, float, float]:
    try:
        records = navigation.loc[satellite]
    except KeyError as error:
        raise ValueError("navigation file has no ephemeris for {}".format(satellite)) from error
    if isinstance(records, pd.Series):
        records = records.to_frame().T

    elapsed_gps = (timestamp - _GPS_EPOCH).total_seconds()
    gps_week = int(elapsed_gps // _GPS_WEEK_SECONDS)
    time_of_week = elapsed_gps - gps_week * _GPS_WEEK_SECONDS
    candidate_rows = []
    for _, record in records.iterrows():
        if float(record["sv_health"]) != 0:
            continue
        ephemeris_week = int(record["gps_week"])
        tk = _wrap_gps_week(
            (gps_week - ephemeris_week) * _GPS_WEEK_SECONDS
            + time_of_week
            - float(record["toe"])
        )
        candidate_rows.append((abs(tk), tk, record))
    if not candidate_rows:
        raise ValueError("navigation file has no healthy ephemeris for {}".format(satellite))

    _, tk, eph = min(candidate_rows, key=lambda item: item[0])
    sqrt_a = float(eph["sqrt_a"])
    eccentricity = float(eph["ecc"])
    semi_major_axis = sqrt_a * sqrt_a
    mean_motion = math.sqrt(_EARTH_GRAVITATIONAL_CONSTANT / semi_major_axis ** 3)
    mean_anomaly = float(eph["m0"]) + (
        mean_motion + float(eph["delta_n"])
    ) * tk
    eccentric_anomaly = mean_anomaly
    for _ in range(12):
        updated = mean_anomaly + eccentricity * math.sin(eccentric_anomaly)
        if abs(updated - eccentric_anomaly) < 1e-13:
            eccentric_anomaly = updated
            break
        eccentric_anomaly = updated

    sin_e = math.sin(eccentric_anomaly)
    cos_e = math.cos(eccentric_anomaly)
    true_anomaly = math.atan2(
        math.sqrt(1 - eccentricity * eccentricity) * sin_e,
        cos_e - eccentricity,
    )
    argument = true_anomaly + float(eph["omega"])
    twice_argument = 2 * argument
    argument += (
        float(eph["cus"]) * math.sin(twice_argument)
        + float(eph["cuc"]) * math.cos(twice_argument)
    )
    radius = (
        semi_major_axis * (1 - eccentricity * cos_e)
        + float(eph["crs"]) * math.sin(twice_argument)
        + float(eph["crc"]) * math.cos(twice_argument)
    )
    inclination = (
        float(eph["inc"])
        + float(eph["idot"]) * tk
        + float(eph["cis"]) * math.sin(twice_argument)
        + float(eph["cic"]) * math.cos(twice_argument)
    )

    orbital_x = radius * math.cos(argument)
    orbital_y = radius * math.sin(argument)
    toe = float(eph["toe"])
    longitude = (
        float(eph["omega0"])
        + (float(eph["dot_omega"]) - _EARTH_ROTATION_RATE) * tk
        - _EARTH_ROTATION_RATE * toe
    )
    x = orbital_x * math.cos(longitude) - orbital_y * math.cos(inclination) * math.sin(longitude)
    y = orbital_x * math.sin(longitude) + orbital_y * math.cos(inclination) * math.cos(longitude)
    z = orbital_y * math.sin(inclination)

    if pseudorange is not None and pseudorange > 0:
        travel_time = pseudorange / _SPEED_OF_LIGHT
        rotation = _EARTH_ROTATION_RATE * travel_time
        x, y = (
            x * math.cos(rotation) + y * math.sin(rotation),
            -x * math.sin(rotation) + y * math.cos(rotation),
        )
    return x, y, z


def _azimuth_elevation(
    satellite_xyz: Tuple[float, float, float],
    receiver_xyz: Tuple[float, float, float],
) -> Tuple[float, float]:
    x, y, z = receiver_xyz
    longitude = math.atan2(y, x)
    latitude = math.atan2(z, math.hypot(x, y) * (1 - 6.69437999014e-3))
    dx, dy, dz = (satellite_xyz[i] - receiver_xyz[i] for i in range(3))
    east = -math.sin(longitude) * dx + math.cos(longitude) * dy
    north = (
        -math.sin(latitude) * math.cos(longitude) * dx
        - math.sin(latitude) * math.sin(longitude) * dy
        + math.cos(latitude) * dz
    )
    up = (
        math.cos(latitude) * math.cos(longitude) * dx
        + math.cos(latitude) * math.sin(longitude) * dy
        + math.sin(latitude) * dz
    )
    azimuth = math.degrees(math.atan2(east, north))
    elevation = math.degrees(math.atan2(up, math.hypot(east, north)))
    return azimuth, elevation


def write_azimuth_elevation_files(
    observation_path: Path,
    navigation_path: Path,
    output_dir: Optional[Path] = None,
    sample_interval: int = 15,
) -> Tuple[Path, Path]:
    """Write TEQC-compatible compact GPS azimuth and elevation files."""
    if sample_interval <= 0:
        raise ValueError("sample interval must be a positive number of seconds")

    output_dir = DATA_2_OUTPUT_DIR if output_dir is None else Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    azimuth_path = output_dir / (observation_path.stem + ".azi")
    elevation_path = output_dir / (observation_path.stem + ".ele")
    _, navigation = read_navigation_dataframe(navigation_path)
    navigation_satellites = sorted(
        satellite
        for satellite, records in navigation.groupby(level="satellite")
        if satellite.startswith("G") and (records["sv_health"] == 0).any()
    )

    with open_observation_file(observation_path) as (header, epochs):
        receiver_xyz = header.approximate_position_xyz
        if receiver_xyz is None:
            raise ValueError("observation header is missing APPROX POSITION XYZ")
        first_epoch = None
        previous_satellites = None
        with (
            azimuth_path.open("w", encoding="ascii", newline="\n") as azimuth_file,
            elevation_path.open("w", encoding="ascii", newline="\n") as elevation_file,
        ):
            for output_file in (azimuth_file, elevation_file):
                output_file.write("COMPACT2\n")
                output_file.write("T_SAMP    {:0.4f}\n".format(sample_interval))
            for epoch in epochs:
                if first_epoch is None:
                    first_epoch = epoch.timestamp
                    mjd = (
                        epoch.timestamp - datetime(1858, 11, 17)
                    ).total_seconds() / 86400.0
                    for output_file in (azimuth_file, elevation_file):
                        output_file.write("START_TIME_MJD  {:0.9f}\n".format(mjd))
                elapsed = (epoch.timestamp - first_epoch).total_seconds()
                if round(elapsed) % sample_interval != 0:
                    continue

                observed_satellites = set(epoch.observations)
                azimuths = []
                elevations = []
                for satellite in navigation_satellites:
                    pseudorange = None
                    if satellite in observed_satellites:
                        for observation_type, observation in zip(
                            header.observation_types,
                            epoch.observations[satellite],
                        ):
                            if observation_type in ("P1", "C1", "P2", "C2"):
                                if observation.value is not None and observation.value > 0:
                                    pseudorange = observation.value
                                    break
                    position = _satellite_position(
                        epoch.timestamp, satellite, navigation, pseudorange
                    )
                    azimuth, elevation = _azimuth_elevation(position, receiver_xyz)
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


def calculate_tec_dataframe(
    observation_path: Path,
    elevation_path: Path,
    shell_height_km: float = 350.0,
) -> pd.DataFrame:
    """Adapt RINEX 2 observations to the reusable TEC calculation API."""
    _, elevations_by_epoch = elevations_from_compact_file(elevation_path)
    elevation_records = [
        {"epoch": timestamp, "satellite": satellite, "elevation_deg": elevation}
        for timestamp, satellites in elevations_by_epoch.items()
        for satellite, elevation in satellites.items()
    ]
    observation_records = []
    with open_observation_file(observation_path) as (header, epochs):
        code_pairs = [
            pair for pair in (("P1", "P2"), ("C1", "C2"))
            if all(code in header.observation_types for code in pair)
        ]
        if not code_pairs:
            raise ValueError(
                "observation file must contain P1/P2 or C1/C2 code measurements"
            )
        code_pair_indices = [
            (pair, tuple(header.observation_types.index(code) for code in pair))
            for pair in code_pairs
        ]
        if "L1" not in header.observation_types or "L2" not in header.observation_types:
            raise ValueError(
                "observation file must contain both L1 and L2 carrier phases "
                "for phase smoothing"
            )
        l1_index = header.observation_types.index("L1")
        l2_index = header.observation_types.index("L2")
        for epoch in epochs:
            for satellite, observations in epoch.observations.items():
                if not satellite.startswith("G"):
                    continue
                selected_pair = code_pairs[0]
                code_1 = code_2 = None
                for pair, (p1_index, p2_index) in code_pair_indices:
                    first = observations[p1_index].value
                    second = observations[p2_index].value
                    if (
                        first is not None
                        and second is not None
                        and first > 0
                        and second > 0
                    ):
                        selected_pair = pair
                        code_1, code_2 = first, second
                        break
                l1_observation = observations[l1_index]
                l2_observation = observations[l2_index]
                observation_records.append({
                    "epoch": epoch.timestamp,
                    "satellite": satellite,
                    "code_1_m": code_1,
                    "code_2_m": code_2,
                    "phase_1_cycles": l1_observation.value,
                    "phase_2_cycles": l2_observation.value,
                    "frequency_1_hz": 1575.42e6,
                    "frequency_2_hz": 1227.60e6,
                    "lli_1": l1_observation.loss_of_lock_indicator,
                    "lli_2": l2_observation.loss_of_lock_indicator,
                    "code_pair": "/".join(selected_pair),
                })

    return calculate_normalized_tec_dataframe(
        pd.DataFrame.from_records(
            observation_records, columns=OBSERVATION_COLUMNS
        ),
        pd.DataFrame.from_records(elevation_records, columns=ELEVATION_COLUMNS),
        shell_height_km=shell_height_km,
        time_system=header.time_system,
    )


def _iter_epochs(stream, header: RinexHeader) -> Iterator[Epoch]:
    pending_line = None
    while True:
        line = pending_line if pending_line is not None else stream.readline()
        pending_line = None
        if not line:
            return
        if not line.strip():
            continue
        if line[60:80].strip() in _HEADER_LABELS:
            continue
        if len(line) < 32:
            if line[:28].strip() == "":
                pending_line = _skip_splice_metadata(stream)
                if pending_line is not None:
                    continue
            raise ValueError("invalid epoch header: {!r}".format(line.rstrip()))

        try:
            year = int(line[1:3])
            year += 1900 if year >= 80 else 2000
            month = int(line[4:6])
            day = int(line[7:9])
            hour = int(line[10:12])
            minute = int(line[13:15])
            seconds = float(line[15:26])
            flag = int(line[28:29])
            satellite_count = int(line[29:32])
            timestamp = datetime(year, month, day, hour, minute) + timedelta(seconds=seconds)
        except (ValueError, OverflowError) as error:
            if line[:28].strip() == "":
                pending_line = _skip_splice_metadata(stream)
                if pending_line is not None:
                    continue
            raise ValueError("invalid epoch header: {!r}".format(line.rstrip())) from error

        if flag in (2, 3, 4, 5):
            for _ in range(satellite_count):
                if not stream.readline():
                    raise ValueError("unexpected end of file in event records")
            continue

        satellites = _read_satellites(line, satellite_count, stream)
        observations_per_line = 5
        lines_per_satellite = (
            len(header.observation_types) + observations_per_line - 1
        ) // observations_per_line
        satellite_observations = {}

        for satellite in satellites:
            fields = []
            for _ in range(lines_per_satellite):
                observation_line = stream.readline()
                if not observation_line:
                    raise ValueError("unexpected end of file in satellite observations")
                observation_line = observation_line.rstrip("\r\n").ljust(80)
                fields.extend(
                    observation_line[start:start + 16]
                    for start in range(0, observations_per_line * 16, 16)
                )
            satellite_observations[satellite] = tuple(
                _read_observation(field)
                for field in fields[:len(header.observation_types)]
            )

        yield Epoch(timestamp, flag, satellite_observations)


@contextmanager
def open_observation_file(path: Path):
    """Open a RINEX 2 file and provide its header and streaming epoch iterator."""
    with path.open("r", encoding="ascii", errors="replace", newline="") as stream:
        header = _read_header(stream)
        yield header, _iter_epochs(stream, header)


def read_observation_dataframe(path: Path):
    """Return the RINEX header and a wide DataFrame indexed by epoch and satellite.

    Each observation type has ``value``, ``lli``, and ``signal_strength``
    subcolumns. Missing observation fields are represented as NaN.
    """
    records = []
    index = []
    fields = ("value", "lli", "signal_strength")

    with open_observation_file(path) as (header, epochs):
        for epoch in epochs:
            for satellite, observations in epoch.observations.items():
                index.append((epoch.timestamp, satellite))
                record = {}
                for observation_type, observation in zip(
                    header.observation_types, observations
                ):
                    record[(observation_type, "value")] = observation.value
                    record[(observation_type, "lli")] = (
                        observation.loss_of_lock_indicator
                    )
                    record[(observation_type, "signal_strength")] = (
                        observation.signal_strength
                    )
                records.append(record)

    columns = pd.MultiIndex.from_product(
        [header.observation_types, fields], names=("observation", "field")
    )
    dataframe = pd.DataFrame(
        records,
        index=pd.MultiIndex.from_tuples(index, names=("epoch", "satellite")),
        columns=columns,
    )
    return header, dataframe


def summarize_observation_dataframe(path: Path, header: RinexHeader, dataframe) -> str:
    epoch_values = dataframe.index.get_level_values("epoch")
    satellites = dataframe.index.get_level_values("satellite").unique()
    value_columns = dataframe.xs("value", axis="columns", level="field")
    observation_count = int(value_columns.notna().to_numpy().sum())
    lines = [
        "File: {}".format(path),
        "RINEX version: {:.2f}".format(header.version),
        "File type: {}".format(header.file_type or "unknown"),
        "Marker: {}".format(header.marker_name or "unknown"),
        "Observation types: {}".format(", ".join(header.observation_types)),
        "Observation epochs: {}".format(epoch_values.unique().size),
        "Satellite-epochs: {}".format(len(dataframe)),
        "Valid observation values: {}".format(observation_count),
        "Satellites: {}".format(", ".join(sorted(satellites)) or "none"),
        "Time range: {} to {}".format(
            epoch_values.min().isoformat(sep=" ") if len(epoch_values) else "none",
            epoch_values.max().isoformat(sep=" ") if len(epoch_values) else "none",
        ),
    ]
    return "\n".join(lines)


def _csv_output_path(source: Path, filename: Optional[Path]) -> Path:
    output_name = filename.name if filename is not None else source.name + ".csv"
    return DATA_2_OUTPUT_DIR / output_name


def main_navigation(argv=None) -> int:
    argument_parser = argparse.ArgumentParser(
        description="Parse and summarize a RINEX 2 navigation file.",
    )
    argument_parser.add_argument(
        "navigation_file",
        type=Path,
        help="RINEX 2 navigation file",
    )
    argument_parser.add_argument(
        "--csv",
        type=Path,
        metavar="OUTPUT",
        help="CSV filename in Data_2/output (default: input filename plus .csv)",
    )
    args = argument_parser.parse_args(argv)

    try:
        header, dataframe = read_navigation_dataframe(args.navigation_file)
        print("File: {}".format(args.navigation_file))
        print("RINEX version: {:.2f}".format(header.version))
        print("File type: {}".format(header.file_type or "unknown"))
        print("Leap seconds: {}".format(header.leap_seconds or "unknown"))
        print("Ephemerides: {}".format(len(dataframe)))
        output_path = _csv_output_path(args.navigation_file, args.csv)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        dataframe.to_csv(output_path)
        print("DataFrame written to: {}".format(output_path))
    except (OSError, ValueError) as error:
        print("Error: {}".format(error), file=sys.stderr)
        return 2
    return 0


def main_angles(argv=None) -> int:
    argument_parser = argparse.ArgumentParser(
        description="Generate compact RINEX 2 azimuth and elevation files.",
    )
    argument_parser.add_argument("observation_file", type=Path)
    argument_parser.add_argument("navigation_file", type=Path)
    argument_parser.add_argument(
        "--interval",
        type=int,
        default=15,
        metavar="SECONDS",
        help="sampling interval in seconds (default: 15)",
    )
    argument_parser.add_argument(
        "--outdir",
        type=Path,
        default=DATA_2_OUTPUT_DIR,
        help="output directory (default: %(default)s)",
    )
    args = argument_parser.parse_args(argv)
    try:
        azimuth_path, elevation_path = write_azimuth_elevation_files(
            args.observation_file,
            args.navigation_file,
            args.outdir,
            args.interval,
        )
        print("Azimuth written to: {}".format(azimuth_path))
        print("Elevation written to: {}".format(elevation_path))
    except (OSError, ValueError) as error:
        print("Error: {}".format(error), file=sys.stderr)
        return 2
    return 0


def main_tec(argv=None) -> int:
    argument_parser = argparse.ArgumentParser(
        description="Calculate uncalibrated GPS code STEC and mapped VTEC.",
    )
    argument_parser.add_argument("observation_file", type=Path)
    argument_parser.add_argument(
        "--ele",
        type=Path,
        help="compact elevation file (default: Data_2/output/<obs-stem>.ele)",
    )
    argument_parser.add_argument(
        "--shell-height",
        type=float,
        default=350.0,
        metavar="KM",
        help="thin-shell ionosphere height in kilometers (default: 350)",
    )
    argument_parser.add_argument(
        "--out",
        type=Path,
        help="output CSV path (default: Data_2/output/<obs-stem>_tec.csv)",
    )
    args = argument_parser.parse_args(argv)
    elevation_path = args.ele or DATA_2_OUTPUT_DIR / (
        args.observation_file.stem + ".ele"
    )
    output_path = args.out or DATA_2_OUTPUT_DIR / (
        args.observation_file.stem + "_tec.csv"
    )
    try:
        dataframe = calculate_tec_dataframe(
            args.observation_file,
            elevation_path,
            args.shell_height,
        )
        if dataframe.empty:
            raise ValueError(
                "no matching epochs with valid dual-frequency code, carrier "
                "phase, and elevation measurements"
            )
        output_path.parent.mkdir(parents=True, exist_ok=True)
        dataframe.to_csv(output_path, index=False)
        print("TEC rows: {}".format(len(dataframe)))
        print("Shell height: {:.1f} km".format(args.shell_height))
        print("Output written to: {}".format(output_path))
        print(
            "Generate plots with: python3 rinex_visualization.py tec {}".format(
                output_path
            )
        )
        print(
            "Uncalibrated STEC range: {:.2f} to {:.2f} TECU".format(
                dataframe["stec_tecu_uncalibrated"].min(),
                dataframe["stec_tecu_uncalibrated"].max(),
            )
        )
        print(
            "Uncalibrated VTEC range: {:.2f} to {:.2f} TECU".format(
                dataframe["phase_smoothed_vtec_tecu_uncalibrated"].min(),
                dataframe["phase_smoothed_vtec_tecu_uncalibrated"].max(),
            )
        )
    except (OSError, ValueError) as error:
        print("Error: {}".format(error), file=sys.stderr)
        return 2
    return 0


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] == "nav":
        return main_navigation(argv[1:])
    if argv and argv[0] == "angles":
        return main_angles(argv[1:])
    if argv and argv[0] == "tec":
        return main_tec(argv[1:])

    default_file = DATA_2_DIR / "kepa3150.13o"
    argument_parser = argparse.ArgumentParser(
        description="Parse and summarize a RINEX 2 observation file."
    )
    argument_parser.add_argument(
        "observation_file",
        nargs="?",
        type=Path,
        default=default_file,
        help="RINEX 2 observation file (default: %(default)s)",
    )
    argument_parser.add_argument(
        "--csv",
        type=Path,
        metavar="OUTPUT",
        help="CSV filename in Data_2/output (default: input filename plus .csv)",
    )
    args = argument_parser.parse_args(argv)

    try:
        header, dataframe = read_observation_dataframe(args.observation_file)
        print(summarize_observation_dataframe(args.observation_file, header, dataframe))
        output_path = _csv_output_path(args.observation_file, args.csv)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        dataframe.to_csv(output_path)
        print("DataFrame written to: {}".format(output_path))
    except (OSError, ValueError) as error:
        print("Error: {}".format(error), file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
