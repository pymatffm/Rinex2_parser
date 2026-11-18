import sys, pathlib
from datetime import date
sys.path.insert(0, str(pathlib.Path(__file__).parents[1]))
sys.path.insert(0, str(pathlib.Path(__file__).parent))
import numpy as np, pandas as pd, pytest
import rinex_3_4_tool as rt, make_sample as ms

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
    assert len(s["plots"]) == 4
    assert all(pathlib.Path(x).stat().st_size > 5000 for x in s["plots"])
    assert len(outputs) == 11
    assert all(path.stem.endswith("_" + stamp) for path in outputs)

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
