import sys, pathlib
from datetime import date
from unittest.mock import patch
sys.path.insert(0, str(pathlib.Path(__file__).parents[1]))
sys.path.insert(0, str(pathlib.Path(__file__).parent))
import numpy as np, pandas as pd, pytest
import rinex_3_4_tool as rt, make_sample as ms
import rinex_visualization as viz
from rinex_visualization import plot_observation_diagnostics


def test_plot_context_uses_station_and_recording_dates():
    observations = pd.DataFrame(
        {
            "time": pd.to_datetime(["2024-03-15 00:00", "2024-03-15 00:02"]),
            "time_utc": pd.to_datetime(["2024-03-14 23:59", "2024-03-15 00:01"]),
        }
    )
    assert viz._observation_plot_context(observations, "", "MERS00TUR") == (
        "MERS00TUR - 2024-03-14 to 2024-03-15"
    )
    assert viz._station_from_filename(
        pathlib.Path("MERS00TUR_R_20240740000_01D_30S_MO.rnx")
    ) == "MERS00TUR"


@pytest.fixture(scope="module", params=[3.04, 4.00])
def obs(request, tmp_path_factory):
    p = tmp_path_factory.mktemp("o") / "t.rnx"
    ms.make_obs(p, ver=request.param)
    return request.param, rt.read_obs(p)

def test_header(obs):
    ver, (h, df) = obs
    assert h.version == ver
    assert h.approx_xyz == pytest.approx(ms.XYZ)
    assert set(h.obs_types) == {"G", "E", "C"} and len(h.obs_types["G"]) == 10
    lat, lon, hgt = h.approx_llh
    assert 49 < lat < 50 and 12 < lon < 13 and 0 < hgt < 1000
    assert h.leap_seconds == 18 and h.time_system == "GPS"

def test_df(obs):
    _, (h, df) = obs
    assert set(df.system) == {"G", "E", "C"}
    assert {"C", "L", "S", "D"} <= set(df.type)
    assert (df[df.type == "S"].value.between(20, 60)).all()
    # pseudorange / phase consistency: C/lambda ≈ L within noise
    w = rt.observations_wide(df, ("C", "L"))
    r = w.dropna(subset=["C1C", "L1C"]); g = r[r.sat.str[0] == "G"]
    assert np.allclose(g.C1C / ms.WL["1"], g.L1C, atol=5)

def test_time_conversion(obs):
    _, (h, df) = obs
    d = (df.time_tai - df.time_gpst).unique(); assert list(d) == [pd.Timedelta(seconds=19)]
    d = (df.time_gpst - df.time_utc).unique(); assert list(d) == [pd.Timedelta(seconds=18)]

def test_scales_roundtrip():
    t = pd.DatetimeIndex(["2016-12-31 23:59:00", "2017-01-01 00:10:00", "2024-03-14"])
    for a in ("UTC", "TAI", "BDT", "GPS"):
        assert (rt.convert_time(rt.convert_time(t, "GPS", a), a, "GPS") == t).all()
    # BDT = GPST - 14 s ; GPST - UTC = 17 s before 2017, 18 s after
    assert (rt.convert_time(t, "GPS", "BDT") - t).total_seconds().tolist() == [-14] * 3
    assert (t - rt.convert_time(t, "GPS", "UTC")).total_seconds().tolist() == [17, 18, 18]
    # BDT 2006-01-01 00:00:00 == GPST 2006-01-01 00:00:14 ; GPS epoch week 0
    wt = rt.gps_week_tow(pd.DatetimeIndex(["1980-01-06", "2024-03-14"]))
    assert wt.week.tolist() == [0, 2305] and wt.tow[1] == 4 * 86400

@pytest.mark.parametrize("ver", [3.04, 4.00])
def test_nav(tmp_path, ver):
    p = tmp_path / "n.rnx"; ms.make_nav(p, ver)
    nav = rt.read_nav(p)
    e = nav.ephemeris
    assert sorted(e.sat) == ["C06", "E11", "G04"]
    assert e.loc[e.sat == "G04", "sqrt_a"].iloc[0] == pytest.approx(5153.8)
    assert e.loc[e.sat == "C06", "toc_gpst"].iloc[0] == pd.Timestamp("2024-03-14 02:00:00") + pd.Timedelta(seconds=14)
    if ver >= 4:
        assert nav.sto.label.iloc[0] == "GAUT" and nav.sto.a0.iloc[0] == pytest.approx(-1.862645149231e-09)
        assert nav.sto.t_tm.iloc[0] == 173500 and nav.sto.a1.iloc[0] == 0
        # Klobuchar ION: alpha0..3 + beta0..3 = 8 values
        assert len(nav.ion.params.iloc[0]) == 8
    else:
        assert nav.header["iono"]["GPSA"][0] == pytest.approx(1.1176e-8)


def test_nav_header_time_correction_with_adjacent_exponents(tmp_path):
    p = tmp_path / "time-correction.rnx"
    p.write_text(
        ms.hdr("     3.04           NAVIGATION DATA     M", "RINEX VERSION / TYPE")
        + ms.hdr(
            "GPGA  3.7834979594D-09-4.440892099D-15  86400 2248",
            "TIME SYSTEM CORR",
        )
        + ms.hdr("", "END OF HEADER")
    )

    correction = rt.read_nav(p).header["time_corr"]["GPGA"]
    assert correction["a0"] == pytest.approx(3.7834979594e-9)
    assert correction["a1"] == pytest.approx(-4.440892099e-15)
    assert correction["t"] == 86400
    assert correction["w"] == 2248


def test_rinex2_rejected(tmp_path):
    p = tmp_path / "o.rnx"
    p.write_text("     2.11           OBSERVATION DATA    G (GPS)             RINEX VERSION / TYPE\n")
    with pytest.raises(ValueError): rt.read_obs(p)

def test_export(tmp_path):
    p = tmp_path / "t.rnx"; ms.make_obs(p, ver=4.00)
    s = rt.export_all(p, tmp_path / "out")
    out = tmp_path / "out" / "v400"
    stamp = date.today().strftime("%d%m%y")
    outputs = list(out.iterdir())
    assert "plots" not in s
    assert len(outputs) == 7
    assert all(path.stem.endswith("_" + stamp) for path in outputs)
    observations = pd.read_csv(
        out / f"observations_long_{stamp}.csv",
        parse_dates=["time", "time_gpst", "time_utc", "time_tai"],
    )
    plots = plot_observation_diagnostics(observations, tmp_path / "plots")
    assert len(plots) == 4
    assert all(path.stat().st_size > 5000 for path in plots)

def test_nav_cli_writes_dated_files_in_version_directory(tmp_path):
    p = tmp_path / "n.rnx"; ms.make_nav(p, 4.00)
    out = tmp_path / "output"
    rt.main(["nav", str(p), "--out", str(out)])
    stamp = date.today().strftime("%d%m%y")
    assert {path.name for path in (out / "v400").iterdir()} == {
        f"ephemeris_{stamp}.csv",
        f"sto_{stamp}.csv",
        f"ion_{stamp}.csv",
    }


@pytest.mark.parametrize("ver", [3.04, 4.00])
def test_generate_compact_azimuth_elevation(tmp_path, ver):
    observation_path = tmp_path / "sample.obs"
    navigation_path = tmp_path / "sample.nav"
    ms.make_obs(observation_path, ver=ver, hours=0.25)
    ms.make_nav(navigation_path, ver=ver)

    receiver = np.asarray(ms.XYZ)
    above_receiver = tuple(receiver * (1 + 20_000_000 / np.linalg.norm(receiver)))
    with patch(
        "rinex_3_4_tool._satellite_position_gnss",
        return_value=above_receiver,
    ):
        azimuth_path, elevation_path = rt.write_azimuth_elevation_files(
            observation_path,
            navigation_path,
            output_dir=tmp_path / "output",
            sample_interval=30,
        )

    assert azimuth_path.is_file()
    assert elevation_path.is_file()
    assert azimuth_path.read_text().startswith("COMPACT2\nT_SAMP    30.0000\n")
    assert elevation_path.read_text().startswith("COMPACT2\nT_SAMP    30.0000\n")
    assert "START_TIME_MJD" in elevation_path.read_text()
    assert len(elevation_path.read_text().splitlines()) > 6
    markers = elevation_path.read_text().splitlines()[3::2]
    assert any(int(marker.split()[0]) > 0 for marker in markers if marker != "-1")


def test_angles_only_processes_requested_observed_constellation(tmp_path):
    observation_path = tmp_path / "sample.obs"
    navigation_path = tmp_path / "sample.nav"
    ms.make_obs(observation_path, ver=4.00, hours=0.25)
    ms.make_nav(navigation_path, ver=4.00)
    receiver = np.asarray(ms.XYZ)
    above_receiver = tuple(receiver * (1 + 20_000_000 / np.linalg.norm(receiver)))

    with patch(
        "rinex_3_4_tool._satellite_position_gnss",
        return_value=above_receiver,
    ) as position:
        rt.write_azimuth_elevation_files(
            observation_path,
            navigation_path,
            output_dir=tmp_path / "output",
            sample_interval=30,
            systems="E",
        )

    assert position.call_count > 0
    assert all(call.args[1].startswith("E") for call in position.call_args_list)
    assert all(len(call.args) == 5 for call in position.call_args_list)


@pytest.mark.parametrize("ver", [3.04, 4.00])
def test_broadcast_orbit_position_is_finite(tmp_path, ver):
    observation_path = tmp_path / "sample.obs"
    navigation_path = tmp_path / "sample.nav"
    ms.make_obs(observation_path, ver=ver, hours=0.25)
    ms.make_nav(navigation_path, ver=ver)
    nav = rt.read_nav(navigation_path)
    timestamp_gpst = pd.Timestamp("2024-03-14 00:00:00")

    position = rt._satellite_position_gnss(
        timestamp_gpst,
        "G04",
        nav.ephemeris,
        23_000_000.0,
    )
    prepared_position = rt._satellite_position_gnss(
        timestamp_gpst,
        "G04",
        nav.ephemeris,
        23_000_000.0,
        rt._prepare_ephemerides(nav.ephemeris)["G04"],
    )

    assert np.isfinite(position).all()
    assert 20_000_000 < np.linalg.norm(position) < 30_000_000
    assert prepared_position == pytest.approx(position)


def test_angles_cli_writes_versioned_outputs(tmp_path):
    observation_path = tmp_path / "sample.obs"
    navigation_path = tmp_path / "sample.nav"
    ms.make_obs(observation_path, ver=4.00, hours=0.25)
    ms.make_nav(navigation_path, ver=4.00)

    assert rt.main([
        "angles",
        str(observation_path),
        str(navigation_path),
        "--out",
        str(tmp_path / "output"),
        "--interval",
        "30",
    ]) == 0
    version_dir = tmp_path / "output" / "v400"
    assert (version_dir / "sample.azi").is_file()
    assert (version_dir / "sample.ele").is_file()


@pytest.mark.parametrize("ver", [3.04, 4.00])
def test_tec_cli_calculates_from_compact_elevation(tmp_path, ver):
    observation_path = tmp_path / "sample.obs"
    navigation_path = tmp_path / "sample.nav"
    ms.make_obs(observation_path, ver=ver, hours=0.25)
    ms.make_nav(navigation_path, ver=ver)
    output_root = tmp_path / "output"
    receiver = np.asarray(ms.XYZ)
    above_receiver = tuple(receiver * (1 + 20_000_000 / np.linalg.norm(receiver)))

    with patch(
        "rinex_3_4_tool._satellite_position_gnss",
        return_value=above_receiver,
    ):
        assert rt.main([
            "angles",
            str(observation_path),
            str(navigation_path),
            "--out",
            str(output_root),
            "--interval",
            "30",
        ]) == 0

    assert rt.main([
        "tec",
        str(observation_path),
        "--out",
        str(output_root),
    ]) == 0
    version_dir = output_root / ("v304" if ver == 3.04 else "v400")
    tec_path = version_dir / "sample_tec.csv"
    tec = pd.read_csv(tec_path, parse_dates=["epoch"])
    assert not tec.empty
    assert {"stec_tecu_uncalibrated", "phase_smoothed_vtec_tecu_uncalibrated"} <= set(
        tec.columns
    )
    assert tec["epoch"].iloc[0] == pd.Timestamp("2024-03-13 23:59:42")
    assert tec["epoch"].iloc[1] == pd.Timestamp("2024-03-14 00:00:12")
