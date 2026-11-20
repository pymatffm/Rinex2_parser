#!/usr/bin/env python3
"""Draw a geographic station map, optionally including ionospheric pierce points.

The map uses a simple longitude/latitude grid and only requires NumPy and
Matplotlib. Optional angle input is a CSV containing ``elevation_deg`` and
``azimuth_deg`` columns, with azimuth clockwise from true north.
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

EARTH_RADIUS_M = 6_378_137.0  # WGS-84 equatorial radius
DEFAULT_SHELL_HEIGHT_KM = 350.0


def floatornan(value) -> float:
    """Convert a nonempty value to float, returning NaN for blank values."""
    if value is None or not str(value).strip():
        return float("nan")
    return float(value)


def digitorzero(value) -> int:
    """Convert a nonempty value to int, returning zero for blank values."""
    if value is None or not str(value).strip():
        return 0
    return int(value)


def padline(line: str, width: int = 16) -> str:
    """Pad a line with spaces to the next multiple of ``width``."""
    if width <= 0:
        raise ValueError("width must be positive")
    padding = (-len(line)) % width
    return line + " " * padding


def ionospheric_pierce_point(
    receiver_latitude_deg: float,
    receiver_longitude_deg: float,
    elevation_deg: float,
    azimuth_deg: float,
    shell_height_km: float = DEFAULT_SHELL_HEIGHT_KM,
) -> Tuple[float, float]:
    """Return thin-shell pierce-point latitude/longitude in degrees.

    Azimuth is measured clockwise from north. Elevation must be above the
    horizon; the shell height is measured above the spherical Earth surface.
    """
    values = (
        receiver_latitude_deg,
        receiver_longitude_deg,
        elevation_deg,
        azimuth_deg,
        shell_height_km,
    )
    if not all(math.isfinite(float(value)) for value in values):
        raise ValueError("coordinates, angles, and shell height must be finite")
    if not -90 <= receiver_latitude_deg <= 90:
        raise ValueError("receiver latitude must be between -90 and 90 degrees")
    if not 0 <= elevation_deg <= 90:
        raise ValueError("satellite elevation must be between 0 and 90 degrees")
    if shell_height_km <= 0:
        raise ValueError("shell height must be positive")

    coefficient = EARTH_RADIUS_M / (EARTH_RADIUS_M + shell_height_km * 1000)
    receiver_latitude = np.deg2rad(receiver_latitude_deg)
    receiver_longitude = np.deg2rad(receiver_longitude_deg)
    elevation = np.deg2rad(elevation_deg)
    azimuth = np.deg2rad(azimuth_deg % 360)
    central_angle = np.pi / 2 - elevation - np.arcsin(
        np.clip(coefficient * np.cos(elevation), -1.0, 1.0)
    )
    pierce_latitude = np.arcsin(
        np.clip(
            np.sin(receiver_latitude) * np.cos(central_angle)
            + np.cos(receiver_latitude) * np.sin(central_angle) * np.cos(azimuth),
            -1.0,
            1.0,
        )
    )

    crosses_pole = (
        receiver_latitude_deg > 70
        and np.tan(central_angle) * np.cos(azimuth)
        > np.tan(np.pi / 2 - receiver_latitude)
    ) or (
        receiver_latitude_deg < -70
        and np.tan(central_angle) * np.cos(azimuth + np.pi)
        > np.tan(np.pi / 2 + receiver_latitude)
    )
    longitude_offset = np.arcsin(
        np.clip(
            np.sin(central_angle) * np.sin(azimuth) / np.cos(pierce_latitude),
            -1.0,
            1.0,
        )
    )
    pierce_longitude = receiver_longitude + longitude_offset
    if crosses_pole:
        pierce_longitude += np.pi

    return (
        float(np.rad2deg(pierce_latitude)),
        float((np.rad2deg(pierce_longitude) + 180) % 360 - 180),
    )


def cartesian_to_ellipsoidal(
    x_m: float, y_m: float, z_m: float
) -> Tuple[float, float, float]:
    """Convert Earth-centered Cartesian coordinates to WGS-84 lon, lat, height."""
    if not all(math.isfinite(value) for value in (x_m, y_m, z_m)):
        raise ValueError("Cartesian coordinates must be finite")
    semi_major = 6_378_137.0
    semi_minor = 6_356_752.314
    eccentricity_squared = 1 - semi_minor**2 / semi_major**2
    longitude = math.atan2(y_m, x_m)
    horizontal_radius = math.hypot(x_m, y_m)
    if horizontal_radius == 0:
        if z_m == 0:
            raise ValueError("the Earth-centered origin has no geodetic coordinates")
        return 0.0, math.copysign(90.0, z_m), abs(z_m) - semi_minor
    latitude = math.atan2(
        z_m, horizontal_radius * (1 - eccentricity_squared)
    )
    for _ in range(10):
        prime_vertical = semi_major / math.sqrt(
            1 - eccentricity_squared * math.sin(latitude) ** 2
        )
        height = horizontal_radius / max(abs(math.cos(latitude)), 1e-15) - prime_vertical
        updated_latitude = math.atan2(
            z_m,
            horizontal_radius
            * (1 - eccentricity_squared * prime_vertical / (prime_vertical + height)),
        )
        if abs(updated_latitude - latitude) < 1e-13:
            latitude = updated_latitude
            break
        latitude = updated_latitude
    prime_vertical = semi_major / math.sqrt(
        1 - eccentricity_squared * math.sin(latitude) ** 2
    )
    height = horizontal_radius / max(abs(math.cos(latitude)), 1e-15) - prime_vertical
    return math.degrees(longitude), math.degrees(latitude), height


def _unwrap_longitudes(
    longitudes: Sequence[float], center_longitude: float
) -> np.ndarray:
    longitudes = np.asarray(longitudes, dtype=float)
    return center_longitude + (longitudes - center_longitude + 180) % 360 - 180


def draw_station_map(
    latitude_deg: float,
    longitude_deg: float,
    output_path: Path = Path("map.png"),
    *,
    latitude_span_deg: float = 30.0,
    longitude_span_deg: float = 60.0,
    pierce_points: Optional[Sequence[Tuple[float, float]]] = None,
    title: Optional[str] = None,
) -> Path:
    """Save a simple geographic map centered on a receiver and return its path."""
    if not math.isfinite(latitude_deg) or not -90 <= latitude_deg <= 90:
        raise ValueError("latitude must be finite and between -90 and 90 degrees")
    if not math.isfinite(longitude_deg):
        raise ValueError("longitude must be finite")
    if not 0 < latitude_span_deg <= 180 or not 0 < longitude_span_deg <= 360:
        raise ValueError("map spans must be positive and no greater than 180/360 degrees")

    center_longitude = (longitude_deg + 180) % 360 - 180
    longitude_min = center_longitude - longitude_span_deg / 2
    longitude_max = center_longitude + longitude_span_deg / 2
    latitude_min = max(-90.0, latitude_deg - latitude_span_deg / 2)
    latitude_max = min(90.0, latitude_deg + latitude_span_deg / 2)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    figure, axis = plt.subplots(figsize=(10, 6), constrained_layout=True)
    axis.set_xlim(longitude_min, longitude_max)
    axis.set_ylim(latitude_min, latitude_max)
    axis.set_xlabel("Longitude (degrees)")
    axis.set_ylabel("Latitude (degrees)")
    axis.set_title(title or "Station and ionospheric pierce points")
    axis.set_facecolor("#d9eff7")
    axis.set_xticks(np.arange(
        math.ceil(longitude_min / 10) * 10,
        longitude_max + 1,
        10,
    ))
    axis.set_yticks(np.arange(
        math.ceil(latitude_min / 10) * 10,
        latitude_max + 1,
        10,
    ))
    axis.grid(color="white", linewidth=1.0)
    axis.scatter(
        center_longitude,
        latitude_deg,
        marker="^",
        s=100,
        color="black",
        label="Receiver",
        zorder=3,
    )

    if pierce_points:
        points = np.asarray(pierce_points, dtype=float)
        if points.ndim != 2 or points.shape[1] != 2:
            raise ValueError("pierce points must be (latitude, longitude) pairs")
        point_lats = points[:, 0]
        point_lons = _unwrap_longitudes(points[:, 1], center_longitude)
        valid = (
            np.isfinite(point_lats)
            & np.isfinite(point_lons)
            & (point_lats >= -90)
            & (point_lats <= 90)
        )
        axis.plot(
            point_lons[valid],
            point_lats[valid],
            color="tab:blue",
            linewidth=1,
            marker=".",
            markersize=3,
            label="Pierce-point track",
        )
    axis.legend(loc="best")
    figure.savefig(output_path, dpi=150)
    plt.close(figure)
    return output_path


def plot_map_globals(
    dphi: float,
    dlambda: float,
    output_path: Path = Path("global_map.png"),
) -> Path:
    """Draw a standalone global longitude/latitude map with a receiver marker."""
    if not math.isfinite(dphi) or not -90 <= dphi <= 90:
        raise ValueError("latitude must be finite and between -90 and 90 degrees")
    if not math.isfinite(dlambda):
        raise ValueError("longitude must be finite")
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure, axis = plt.subplots(figsize=(12, 6), constrained_layout=True)
    axis.set_xlim(-180, 180)
    axis.set_ylim(-90, 90)
    axis.set_xticks(np.arange(-180, 181, 30))
    axis.set_yticks(np.arange(-90, 91, 15))
    axis.set_facecolor("#d9eff7")
    axis.grid(color="white", linewidth=0.8)
    axis.set_xlabel("Longitude (degrees)")
    axis.set_ylabel("Latitude (degrees)")
    axis.set_title("Global receiver location")
    axis.scatter(
        (dlambda + 180) % 360 - 180,
        dphi,
        marker="^",
        s=100,
        color="black",
        label="Receiver",
        zorder=3,
    )
    axis.legend(loc="best")
    figure.savefig(output_path, dpi=150)
    plt.close(figure)
    return output_path


def funky_earth_plot(
    dphi: float,
    dlambda: float,
    output_path: Path = Path("funky_earth.png"),
) -> Path:
    """Draw a standalone globe-style orthographic view centered on the receiver."""
    if not math.isfinite(dphi) or not -90 <= dphi <= 90:
        raise ValueError("latitude must be finite and between -90 and 90 degrees")
    if not math.isfinite(dlambda):
        raise ValueError("longitude must be finite")
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    center_latitude = np.deg2rad(dphi)
    center_longitude = np.deg2rad(dlambda)

    def project(latitude, longitude):
        latitude = np.deg2rad(np.asarray(latitude, dtype=float))
        longitude = np.deg2rad(np.asarray(longitude, dtype=float))
        delta_longitude = longitude - center_longitude
        visible = (
            np.sin(center_latitude) * np.sin(latitude)
            + np.cos(center_latitude) * np.cos(latitude) * np.cos(delta_longitude)
        ) >= 0
        x = np.cos(latitude) * np.sin(delta_longitude)
        y = (
            np.cos(center_latitude) * np.sin(latitude)
            - np.sin(center_latitude) * np.cos(latitude) * np.cos(delta_longitude)
        )
        return x, y, visible

    figure, axis = plt.subplots(figsize=(8, 8), constrained_layout=True)
    globe = plt.Circle((0, 0), 1, color="#d9eff7", zorder=0)
    axis.add_patch(globe)
    longitudes = np.linspace(-180, 180, 361)
    latitudes = np.linspace(-90, 90, 181)
    for latitude in np.arange(-75, 90, 15):
        x, y, visible = project(np.full_like(longitudes, latitude), longitudes)
        x[~visible] = np.nan
        y[~visible] = np.nan
        axis.plot(x, y, color="white", linewidth=0.7)
    for longitude in np.arange(-180, 180, 30):
        x, y, visible = project(latitudes, np.full_like(latitudes, longitude))
        x[~visible] = np.nan
        y[~visible] = np.nan
        axis.plot(x, y, color="white", linewidth=0.7)
    axis.scatter(0, 0, marker="^", s=100, color="black", label="Receiver", zorder=3)
    axis.set_aspect("equal")
    axis.set_xlim(-1.05, 1.05)
    axis.set_ylim(-1.05, 1.05)
    axis.axis("off")
    axis.set_title("Receiver-centered globe")
    axis.legend(loc="lower center")
    figure.savefig(output_path, dpi=150)
    plt.close(figure)
    return output_path


def plot_map(
    lat_array,
    lon_array,
    interesting_satellites,
    dphi: float,
    dlambda: float,
    output_path: Path = Path("ipp_map.png"),
    event_index: Optional[int] = None,
) -> Path:
    """Plot satellite pierce-point tracks; arrays are epoch-by-satellite."""
    latitudes = np.asarray(lat_array, dtype=float)
    longitudes = np.asarray(lon_array, dtype=float)
    if latitudes.ndim != 2 or longitudes.shape != latitudes.shape:
        raise ValueError("latitude and longitude arrays must have matching 2-D shapes")
    if len(interesting_satellites) != latitudes.shape[1]:
        raise ValueError("satellite labels must match the track array columns")
    if event_index is not None and not 0 <= event_index < latitudes.shape[0]:
        raise ValueError("event index is outside the track array")
    if not math.isfinite(dphi) or not -90 <= dphi <= 90:
        raise ValueError("latitude must be finite and between -90 and 90 degrees")
    if not math.isfinite(dlambda):
        raise ValueError("longitude must be finite")

    center_longitude = (dlambda + 180) % 360 - 180
    valid_lats = latitudes[
        np.isfinite(latitudes) & (latitudes >= -90) & (latitudes <= 90)
    ]
    valid_lons = _unwrap_longitudes(
        longitudes[np.isfinite(longitudes)], center_longitude
    )
    latitude_values = np.append(valid_lats, dphi)
    longitude_values = np.append(valid_lons, center_longitude)
    latitude_min = max(-90.0, float(latitude_values.min()) - 5)
    latitude_max = min(90.0, float(latitude_values.max()) + 5)
    longitude_min = float(longitude_values.min()) - 5
    longitude_max = float(longitude_values.max()) + 5
    if latitude_max - latitude_min < 10:
        latitude_min = max(-90.0, dphi - 5)
        latitude_max = min(90.0, dphi + 5)
    if longitude_max - longitude_min < 10:
        longitude_min = center_longitude - 5
        longitude_max = center_longitude + 5
    figure, axis = plt.subplots(figsize=(10, 6), constrained_layout=True)
    axis.set_xlim(longitude_min, longitude_max)
    axis.set_ylim(latitude_min, latitude_max)
    axis.set_facecolor("#d9eff7")
    axis.grid(color="white", linewidth=0.8)
    axis.set_xlabel("Longitude (degrees)")
    axis.set_ylabel("Latitude (degrees)")
    axis.set_title("Satellite ionospheric pierce-point tracks")
    axis.scatter(
        center_longitude, dphi, marker="^", s=100, color="black",
        label="Receiver", zorder=3,
    )
    for index, satellite in enumerate(interesting_satellites):
        track_lons = _unwrap_longitudes(longitudes[:, index], center_longitude)
        track_lats = latitudes[:, index]
        valid = (
            np.isfinite(track_lats)
            & np.isfinite(track_lons)
            & (track_lats >= -90)
            & (track_lats <= 90)
        )
        if not valid.any():
            continue
        label = (
            "G{:02d}".format(satellite + 1)
            if isinstance(satellite, (int, np.integer))
            else str(satellite)
        )
        axis.plot(track_lons[valid], track_lats[valid], linewidth=1, label=label)
        if event_index is not None and valid[event_index]:
            axis.scatter(
                track_lons[event_index], track_lats[event_index],
                marker="*", color="red", s=70, zorder=4,
            )
    axis.legend(loc="best", fontsize="small")
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=150)
    plt.close(figure)
    return output_path


def load_angle_csv(
    path: Path,
    latitude_deg: float,
    longitude_deg: float,
    shell_height_km: float = DEFAULT_SHELL_HEIGHT_KM,
) -> list[Tuple[float, float]]:
    """Read an angle CSV and convert its elevation/azimuth rows to pierce points."""
    angles = pd.read_csv(path)
    required = {"elevation_deg", "azimuth_deg"}
    missing = required - set(angles.columns)
    if missing:
        raise ValueError(
            "angle CSV is missing columns: {}".format(", ".join(sorted(missing)))
        )
    points = []
    for elevation, azimuth in angles.loc[
        :, ["elevation_deg", "azimuth_deg"]
    ].itertuples(index=False, name=None):
        if pd.isna(elevation) or pd.isna(azimuth) or elevation < 0:
            continue
        points.append(
            ionospheric_pierce_point(
                latitude_deg,
                longitude_deg,
                float(elevation),
                float(azimuth),
                shell_height_km,
            )
        )
    return points


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--mode",
        choices=("station", "global", "funky-earth", "tracks"),
        default="station",
        help="map style; tracks requires --tracks CSV",
    )
    parser.add_argument("--latitude", type=float, required=True, help="receiver latitude")
    parser.add_argument("--longitude", type=float, required=True, help="receiver longitude")
    parser.add_argument("--output", type=Path, default=Path("map.png"))
    parser.add_argument("--angles", type=Path, help="CSV with elevation_deg and azimuth_deg")
    parser.add_argument(
        "--tracks",
        type=Path,
        help="CSV with satellite, latitude_deg, longitude_deg columns",
    )
    parser.add_argument("--shell-height", type=float, default=DEFAULT_SHELL_HEIGHT_KM)
    parser.add_argument("--latitude-span", type=float, default=30.0)
    parser.add_argument("--longitude-span", type=float, default=60.0)
    parser.add_argument("--title")
    args = parser.parse_args(argv)

    try:
        if args.mode == "global":
            output_path = plot_map_globals(
                args.latitude, args.longitude, args.output
            )
        elif args.mode == "funky-earth":
            output_path = funky_earth_plot(
                args.latitude, args.longitude, args.output
            )
        elif args.mode == "tracks":
            if not args.tracks:
                raise ValueError("--tracks is required with --mode tracks")
            tracks = pd.read_csv(args.tracks)
            required = {"satellite", "latitude_deg", "longitude_deg"}
            missing = required - set(tracks.columns)
            if missing:
                raise ValueError(
                    "track CSV is missing columns: {}".format(
                        ", ".join(sorted(missing))
                    )
                )
            pivot_lat = tracks.pivot(
                columns="satellite", values="latitude_deg"
            )
            pivot_lon = tracks.pivot(
                columns="satellite", values="longitude_deg"
            )
            output_path = plot_map(
                pivot_lat.to_numpy(),
                pivot_lon.to_numpy(),
                list(pivot_lat.columns),
                args.latitude,
                args.longitude,
                args.output,
            )
        else:
            points = (
                load_angle_csv(
                    args.angles,
                    args.latitude,
                    args.longitude,
                    args.shell_height,
                )
                if args.angles
                else None
            )
            output_path = draw_station_map(
                args.latitude,
                args.longitude,
                args.output,
                latitude_span_deg=args.latitude_span,
                longitude_span_deg=args.longitude_span,
                pierce_points=points,
                title=args.title,
            )
    except (OSError, ValueError, pd.errors.ParserError) as error:
        parser.error(str(error))
    print("Map written to: {}".format(output_path))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
