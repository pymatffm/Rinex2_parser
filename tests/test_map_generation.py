import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parents[1]))

import pytest

import map_generation as maps


def test_pierce_point_is_ahead_of_receiver():
    latitude, longitude = maps.ionospheric_pierce_point(
        0.0, 0.0, elevation_deg=45.0, azimuth_deg=90.0
    )

    assert latitude == pytest.approx(0.0)
    assert longitude > 0


def test_cartesian_to_ellipsoidal_round_trip_location():
    longitude, latitude, height = maps.cartesian_to_ellipsoidal(
        4_239_146.6414, 2_886_967.1245, 3_778_874.4800
    )

    assert longitude == pytest.approx(34.3, abs=0.1)
    assert latitude == pytest.approx(36.5664, abs=0.001)
    assert height == pytest.approx(34.7, abs=0.1)


def test_draw_station_map_without_satellite_arrays(tmp_path):
    output = maps.draw_station_map(
        36.2,
        34.3,
        tmp_path / "station.png",
        pierce_points=[(37.0, 35.0), (38.0, 36.0)],
    )

    assert output.is_file()
    assert output.stat().st_size > 1000


@pytest.mark.parametrize(
    "plotter,filename",
    [
        (maps.plot_map_globals, "global.png"),
        (maps.funky_earth_plot, "globe.png"),
    ],
)
def test_standalone_global_plotters(tmp_path, plotter, filename):
    output = plotter(36.2, 34.3, tmp_path / filename)

    assert output.is_file()
    assert output.stat().st_size > 1000


def test_standalone_satellite_track_plotter(tmp_path):
    output = maps.plot_map(
        [[37.0, 38.0], [38.0, 39.0]],
        [[35.0, 36.0], [36.0, 37.0]],
        [0, "E02"],
        36.0,
        34.0,
        tmp_path / "tracks.png",
        event_index=1,
    )

    assert output.is_file()
    assert output.stat().st_size > 1000


def test_standalone_cli_draws_station_map(tmp_path):
    output = tmp_path / "standalone.png"

    assert maps.main([
        "--latitude",
        "36.2",
        "--longitude",
        "34.3",
        "--output",
        str(output),
    ]) == 0
    assert output.is_file()


def test_global_cli_mode(tmp_path):
    output = tmp_path / "global_cli.png"

    assert maps.main([
        "--mode",
        "global",
        "--latitude",
        "36.2",
        "--longitude",
        "34.3",
        "--output",
        str(output),
    ]) == 0
    assert output.is_file()


def test_angle_csv_converts_rows_to_pierce_points(tmp_path):
    angle_file = tmp_path / "angles.csv"
    angle_file.write_text(
        "elevation_deg,azimuth_deg\n45,90\n-1,0\n",
        encoding="utf-8",
    )

    points = maps.load_angle_csv(angle_file, 0.0, 0.0)

    assert len(points) == 1
    assert points[0][1] > 0
