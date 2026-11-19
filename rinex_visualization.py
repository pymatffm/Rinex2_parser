"""Shared plotting and map generation for normalized RINEX parser outputs."""
from __future__ import annotations

import argparse
import sys
from datetime import date, datetime, timedelta
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

SYSTEMS = {
    "G": "GPS",
    "E": "Galileo",
    "C": "BeiDou",
    "R": "GLONASS",
    "J": "QZSS",
    "I": "NavIC",
    "S": "SBAS",
}
_BAND_LABELS = {
    "G": {"1": "L1", "2": "L2", "5": "L5"},
    "E": {"1": "E1", "5": "E5a", "7": "E5b", "8": "E5", "6": "E6"},
    "C": {"1": "B1C", "2": "B1I", "5": "B2a", "7": "B2b", "8": "B2", "6": "B3"},
}
_SYSTEM_COLORS = {
    "G": "#1f77b4",
    "E": "#2ca02c",
    "C": "#d62728",
    "R": "#9467bd",
    "J": "#8c564b",
    "I": "#e377c2",
    "S": "#7f7f7f",
}
_EARTH_RADIUS_KM = 6371.0


def _matplotlib():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def _dated_path(output_dir: Path, filename: str, date_stamp: str) -> Path:
    path = Path(filename)
    return output_dir / "{}_{}{}".format(path.stem, date_stamp, path.suffix)


def plot_vtec_by_satellite(
    dataframe: pd.DataFrame,
    output_dir: Path,
    polynomial_degree: int = 3,
) -> Tuple[Path, ...]:
    """Plot phase-smoothed VTEC and a polynomial fit with residuals per satellite."""
    if polynomial_degree < 0:
        raise ValueError("polynomial degree cannot be negative")
    required = {
        "epoch",
        "satellite",
        "phase_smoothed_vtec_tecu_uncalibrated",
    }
    if not required.issubset(dataframe.columns):
        raise ValueError(
            "TEC CSV must contain {}".format(", ".join(sorted(required)))
        )
    if dataframe.empty:
        raise ValueError("cannot plot VTEC because there are no TEC records")

    plt = _matplotlib()
    from matplotlib.ticker import FormatStrFormatter

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    plot_paths = []
    for satellite, satellite_data in dataframe.groupby("satellite", sort=True):
        satellite_data = satellite_data.sort_values("epoch")
        times = pd.to_datetime(satellite_data["epoch"])
        utc_hours = np.array([
            time.hour + time.minute / 60 + time.second / 3600
            for time in times
        ])
        values = satellite_data[
            "phase_smoothed_vtec_tecu_uncalibrated"
        ].to_numpy(dtype=float)
        if not np.isfinite(values).all():
            raise ValueError("VTEC values for {} must be finite".format(satellite))
        degree = min(polynomial_degree, len(values) - 1)
        coefficients = np.polyfit(utc_hours, values, degree)
        fitted = np.polyval(coefficients, utc_hours)

        figure, (axis, residual_axis) = plt.subplots(
            2,
            1,
            figsize=(10, 7),
            sharex=True,
            gridspec_kw={"height_ratios": [3, 1]},
        )
        axis.plot(utc_hours, values, label="Phase smoothed", linewidth=1.2)
        axis.plot(
            utc_hours,
            fitted,
            label="{}-degree polyfit".format(degree),
            linewidth=1.2,
        )
        axis.set_title("Ionospheric Delay for {}".format(satellite))
        axis.set_ylabel("VTEC (TECU)")
        axis.grid(True, alpha=0.3)
        axis.legend()

        residual_axis.plot(utc_hours, values - fitted, linewidth=0.9)
        residual_axis.axhline(0, color="black", linewidth=0.8, linestyle="--")
        residual_axis.set_ylabel("Residual\n(TECU)")
        residual_axis.set_xlabel("UTC Time (hh.hh)")
        residual_axis.xaxis.set_major_formatter(FormatStrFormatter("%.2f"))
        residual_axis.grid(True, alpha=0.3)

        figure.tight_layout()
        plot_path = output_dir / "{}_vtec.png".format(satellite)
        figure.savefig(plot_path, dpi=150)
        plt.close(figure)
        plot_paths.append(plot_path)

    return tuple(plot_paths)


def _band_label(system: str, band: str) -> str:
    return _BAND_LABELS.get(system, {}).get(band, "{}{}".format(system, band))


def plot_observation_diagnostics(
    dataframe: pd.DataFrame,
    output_dir: Path,
    title: str = "",
    date_stamp: Optional[str] = None,
) -> List[Path]:
    """Create observation diagnostics from the RINEX 3/4 long CSV schema."""
    required = {"time", "time_gpst", "sat", "system", "type", "band", "value"}
    if not required.issubset(dataframe.columns):
        raise ValueError(
            "observation CSV must contain {}".format(", ".join(sorted(required)))
        )
    if dataframe.empty:
        raise ValueError("cannot plot diagnostics because there are no observations")

    plt = _matplotlib()
    import matplotlib.dates as mdates

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    date_stamp = date_stamp or date.today().strftime("%d%m%y")
    colors = _SYSTEM_COLORS
    files = []
    time_formatter = mdates.DateFormatter("%H:%M")

    snr = dataframe[dataframe["type"] == "S"].copy()
    if not snr.empty:
        snr["signal"] = [
            _band_label(system, band)
            for system, band in zip(snr["system"], snr["band"])
        ]
        figure, axis = plt.subplots(figsize=(11, 5))
        line_styles = ("-", "--", ":", "-.")
        signals = sorted(snr["signal"].unique())
        line_style = {signal: line_styles[i % len(line_styles)]
                      for i, signal in enumerate(signals)}
        for (system, signal), group in snr.groupby(["system", "signal"]):
            means = group.groupby("time_gpst")["value"].mean()
            axis.plot(
                means.index,
                means.values,
                linewidth=1.2,
                label="{} {}".format(SYSTEMS.get(system, system), signal),
                color=colors.get(system),
                linestyle=line_style[signal],
            )
        axis.set(
            xlabel="GPST",
            ylabel="Mean C/N0 (dB-Hz)",
            title="Signal-to-noise ratio over time {}".format(title).strip(),
        )
        axis.xaxis.set_major_formatter(time_formatter)
        axis.grid(alpha=0.3)
        axis.legend(ncol=3, fontsize=8)
        path = _dated_path(output_dir, "snr_vs_time.png", date_stamp)
        figure.tight_layout()
        figure.savefig(path, dpi=150)
        plt.close(figure)
        files.append(path)

        figure, axis = plt.subplots(figsize=(11, 5))
        groups = [
            ("{}\n{}".format(system, signal), group["value"].to_numpy())
            for (system, signal), group in snr.groupby(["system", "signal"])
        ]
        axis.boxplot(
            [values for _, values in groups],
            tick_labels=[label for label, _ in groups],
            showfliers=False,
        )
        axis.set(ylabel="C/N0 (dB-Hz)", title="SNR distribution by signal")
        axis.grid(alpha=0.3, axis="y")
        path = _dated_path(output_dir, "snr_distribution.png", date_stamp)
        figure.tight_layout()
        figure.savefig(path, dpi=150)
        plt.close(figure)
        files.append(path)

    tracked = dataframe.drop_duplicates(["time", "sat"])
    visibility = tracked.groupby(["time", "system"])["sat"].nunique().unstack(
        fill_value=0
    )
    visibility["total"] = visibility.sum(axis=1)
    figure, axis = plt.subplots(figsize=(11, 5))
    systems = [system for system in visibility.columns if system != "total"]
    axis.stackplot(
        visibility.index,
        [visibility[system] for system in systems],
        labels=[SYSTEMS.get(system, system) for system in systems],
        colors=[colors.get(system) for system in systems],
        alpha=0.85,
    )
    axis.plot(visibility.index, visibility["total"], "k-", linewidth=0.8, label="Total")
    axis.set(
        xlabel="Observation time (file scale)",
        ylabel="Satellites tracked",
        title="Satellite visibility over time {}".format(title).strip(),
        ylim=(0, max(1, visibility["total"].max() * 1.15)),
    )
    axis.xaxis.set_major_formatter(time_formatter)
    axis.grid(alpha=0.3)
    axis.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.12),
        ncol=4,
        fontsize=8,
        frameon=False,
    )
    path = _dated_path(output_dir, "satellite_count.png", date_stamp)
    figure.tight_layout()
    figure.savefig(path, dpi=150)
    plt.close(figure)
    files.append(path)

    satellites = sorted(tracked["sat"].unique())
    satellite_index = {satellite: index for index, satellite in enumerate(satellites)}
    figure, axis = plt.subplots(
        figsize=(11, max(4, 0.22 * len(satellites) + 1.5))
    )
    axis.scatter(
        tracked["time"],
        tracked["sat"].map(satellite_index),
        s=6,
        marker="|",
        c=[colors.get(system) for system in tracked["system"]],
    )
    axis.set_yticks(range(len(satellites)))
    axis.set_yticklabels(satellites, fontsize=7)
    axis.set(xlabel="Observation time (file scale)", title="Satellite tracking timeline")
    axis.xaxis.set_major_formatter(time_formatter)
    axis.grid(alpha=0.3, axis="x")
    path = _dated_path(output_dir, "satellite_timeline.png", date_stamp)
    figure.tight_layout()
    figure.savefig(path, dpi=150)
    plt.close(figure)
    files.append(path)
    return files


def _read_compact_angles(path: Path):
    with Path(path).open("r", encoding="ascii") as stream:
        if stream.readline().strip() != "COMPACT2":
            raise ValueError("{} is not a COMPACT2 angle file".format(path))
        interval_parts = stream.readline().split()
        start_parts = stream.readline().split()
        if (
            len(interval_parts) != 2
            or interval_parts[0] != "T_SAMP"
            or len(start_parts) != 2
            or start_parts[0] != "START_TIME_MJD"
        ):
            raise ValueError("{} has invalid compact angle metadata".format(path))
        try:
            interval = float(interval_parts[1])
            start_mjd = float(start_parts[1])
        except ValueError as error:
            raise ValueError("{} has invalid compact angle metadata".format(path)) from error
        if interval <= 0:
            raise ValueError("{} has an invalid sampling interval".format(path))

        start_time = datetime(1858, 11, 17) + timedelta(days=start_mjd)
        start_time = (start_time + timedelta(microseconds=500000)).replace(
            microsecond=0
        )
        satellites = None
        epochs = []
        while True:
            marker_line = stream.readline()
            if not marker_line:
                break
            values_line = stream.readline()
            if not values_line:
                raise ValueError("{} has an incomplete angle row".format(path))
            marker = marker_line.split()
            if not marker:
                raise ValueError("{} has an empty satellite marker".format(path))
            if marker[0] != "-1":
                try:
                    count = int(marker[0])
                except ValueError as error:
                    raise ValueError("{} has an invalid satellite count".format(path)) from error
                satellites = marker[1:]
                if len(satellites) != count:
                    raise ValueError("{} satellite count does not match".format(path))
            elif satellites is None:
                raise ValueError("{} uses -1 before declaring satellites".format(path))
            try:
                values = [float(value) for value in values_line.split()]
            except ValueError as error:
                raise ValueError("{} has an invalid angle value".format(path)) from error
            if len(values) != len(satellites):
                raise ValueError("{} angle count does not match satellites".format(path))
            epochs.append(dict(zip(satellites, values)))

    return interval, start_time, epochs


def _pierce_point(
    receiver_latitude: float,
    receiver_longitude: float,
    azimuth_degrees: float,
    elevation_degrees: float,
    shell_height_km: float,
) -> Tuple[float, float]:
    receiver_lat = np.radians(receiver_latitude)
    receiver_lon = np.radians(receiver_longitude)
    azimuth = np.radians(azimuth_degrees % 360)
    elevation = np.radians(elevation_degrees)
    ratio = _EARTH_RADIUS_KM / (_EARTH_RADIUS_KM + shell_height_km)
    central_angle = np.pi / 2 - elevation - np.arcsin(
        np.clip(ratio * np.cos(elevation), -1.0, 1.0)
    )
    latitude = np.arcsin(
        np.sin(receiver_lat) * np.cos(central_angle)
        + np.cos(receiver_lat) * np.sin(central_angle) * np.cos(azimuth)
    )
    longitude_delta = np.arctan2(
        np.sin(central_angle) * np.sin(azimuth) * np.cos(receiver_lat),
        np.cos(central_angle) - np.sin(receiver_lat) * np.sin(latitude),
    )
    longitude = (receiver_lon + longitude_delta + np.pi) % (2 * np.pi) - np.pi
    return float(np.degrees(latitude)), float(np.degrees(longitude))


def plot_pierce_point_map(
    azimuth_path: Path,
    elevation_path: Path,
    receiver_latitude: float,
    receiver_longitude: float,
    output_path: Optional[Path] = None,
    shell_height_km: float = 350.0,
) -> Path:
    """Plot satellite ionospheric pierce points from matching COMPACT2 angle files."""
    if not -90 <= receiver_latitude <= 90:
        raise ValueError("receiver latitude must be between -90 and 90 degrees")
    if not -180 <= receiver_longitude <= 180:
        raise ValueError("receiver longitude must be between -180 and 180 degrees")
    if shell_height_km <= 0:
        raise ValueError("shell height must be positive")

    azimuth_interval, azimuth_start, azimuth_epochs = _read_compact_angles(
        azimuth_path
    )
    elevation_interval, elevation_start, elevation_epochs = _read_compact_angles(
        elevation_path
    )
    if (
        azimuth_interval != elevation_interval
        or azimuth_start != elevation_start
        or len(azimuth_epochs) != len(elevation_epochs)
    ):
        raise ValueError("azimuth and elevation files do not describe the same epochs")

    points: Dict[str, List[Tuple[float, float]]] = {}
    for azimuths, elevations in zip(azimuth_epochs, elevation_epochs):
        if azimuths.keys() != elevations.keys():
            raise ValueError("azimuth and elevation satellite lists do not match")
        for satellite, elevation in elevations.items():
            latitude, longitude = _pierce_point(
                receiver_latitude,
                receiver_longitude,
                azimuths[satellite],
                elevation,
                shell_height_km,
            )
            points.setdefault(satellite, []).append((longitude, latitude))
    if not points:
        raise ValueError("no satellite geometry is available to map")

    plt = _matplotlib()
    from matplotlib.lines import Line2D

    output_path = Path(output_path) if output_path else Path(azimuth_path).with_name(
        Path(azimuth_path).stem + "_pierce_points.png"
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure, axis = plt.subplots(figsize=(12, 6))
    colors = plt.get_cmap("tab20")
    legend = []
    for index, (satellite, coordinates) in enumerate(sorted(points.items())):
        longitude, latitude = np.asarray(coordinates).T
        color = colors(index % 20)
        axis.plot(longitude, latitude, ".", markersize=2, color=color, alpha=0.7)
        legend.append(Line2D([], [], color=color, marker=".", linestyle="",
                             label=satellite))
    axis.scatter(
        [receiver_longitude],
        [receiver_latitude],
        marker="^",
        color="black",
        label="Receiver",
        zorder=5,
    )
    axis.set(
        xlim=(-180, 180),
        ylim=(-90, 90),
        xlabel="Longitude (degrees)",
        ylabel="Latitude (degrees)",
        title="Ionospheric pierce points ({:.0f} km shell)".format(shell_height_km),
    )
    axis.set_xticks(np.arange(-180, 181, 30))
    axis.set_yticks(np.arange(-90, 91, 15))
    axis.grid(alpha=0.3)
    axis.legend(
        handles=legend + [
            Line2D([], [], color="black", marker="^", linestyle="", label="Receiver")
        ],
        loc="upper center",
        bbox_to_anchor=(0.5, -0.12),
        ncol=min(8, len(legend) + 1),
        fontsize=8,
    )
    figure.tight_layout()
    figure.savefig(output_path, dpi=150)
    plt.close(figure)
    return output_path


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description="Generate plots and maps from RINEX parser outputs."
    )
    commands = parser.add_subparsers(dest="command", required=True)

    tec = commands.add_parser("tec", help="plot phase-smoothed VTEC CSV")
    tec.add_argument("csv", type=Path)
    tec.add_argument("--outdir", type=Path)
    tec.add_argument("--degree", type=int, default=3)

    observations = commands.add_parser(
        "observations", help="plot observation diagnostics from RINEX 3/4 long CSV"
    )
    observations.add_argument("csv", type=Path)
    observations.add_argument("--outdir", type=Path)

    map_parser = commands.add_parser(
        "map", help="map ionospheric pierce points from compact azimuth/elevation"
    )
    map_parser.add_argument("--azimuth", type=Path, required=True)
    map_parser.add_argument("--elevation", type=Path, required=True)
    map_parser.add_argument("--latitude", type=float, required=True)
    map_parser.add_argument("--longitude", type=float, required=True)
    map_parser.add_argument("--shell-height", type=float, default=350.0)
    map_parser.add_argument("--out", type=Path)

    args = parser.parse_args(argv)
    try:
        if args.command == "tec":
            dataframe = pd.read_csv(args.csv, parse_dates=["epoch"])
            output_dir = args.outdir or args.csv.with_name(
                args.csv.stem.replace("_tec", "") + "_vtec_plots"
            )
            paths = plot_vtec_by_satellite(
                dataframe, output_dir, polynomial_degree=args.degree
            )
        elif args.command == "observations":
            dataframe = pd.read_csv(
                args.csv, parse_dates=["time", "time_gpst", "time_utc", "time_tai"]
            )
            output_dir = args.outdir or args.csv.parent / "plots"
            paths = plot_observation_diagnostics(dataframe, output_dir)
        else:
            paths = [
                plot_pierce_point_map(
                    args.azimuth,
                    args.elevation,
                    args.latitude,
                    args.longitude,
                    output_path=args.out,
                    shell_height_km=args.shell_height,
                )
            ]
        print("Generated {} visualization file(s):".format(len(paths)))
        for path in paths:
            print(path)
    except (OSError, ValueError, KeyError) as error:
        print("Error: {}".format(error), file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
