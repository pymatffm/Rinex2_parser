#!/usr/bin/env python3
"""Parse and summarize a RINEX 2 observation file."""

import argparse
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime, timedelta
import math
from pathlib import Path
import sys
from typing import Dict, Iterator, Optional, Tuple

import pandas as pd


@dataclass(frozen=True)
class RinexHeader:
    version: float
    file_type: str
    marker_name: str
    observation_types: Tuple[str, ...]


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
    observation_types = []
    observation_type_count = None

    while True:
        line = stream.readline()
        if not line:
            raise ValueError("missing END OF HEADER")
        label = line[60:80].strip()

        if label == "MARKER NAME":
            marker_name = line[:60].strip()
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


def _iter_epochs(stream, header: RinexHeader) -> Iterator[Epoch]:
    pending_line = None
    while True:
        line = pending_line if pending_line is not None else stream.readline()
        pending_line = None
        if not line:
            return
        if not line.strip():
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


def summarize_observation_file(path: Path) -> str:
    header, dataframe = read_observation_dataframe(path)
    return summarize_observation_dataframe(path, header, dataframe)


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


def main(argv=None) -> int:
    default_file = Path(__file__).resolve().parent / "Data" / "kepa3150.13o"
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
        help="write the epoch/satellite observation table to a CSV file",
    )
    args = argument_parser.parse_args(argv)

    try:
        header, dataframe = read_observation_dataframe(args.observation_file)
        print(summarize_observation_dataframe(args.observation_file, header, dataframe))
        if args.csv:
            dataframe.to_csv(args.csv)
            print("DataFrame written to: {}".format(args.csv))
    except (OSError, ValueError) as error:
        print("Error: {}".format(error), file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
