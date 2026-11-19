# Introduction
This repository provides scripts for parsing RINEX 2, RINEX 3, and RINEX 4 files and exporting observation and
navigation data to CSV for all GNSS constellations.

For RINEX 2 data, the parsers can also generate per-satellite elevation and
azimuth files, similar to those produced by the legacy teqc utility. A separate
visualization command generates satellite plots and pierce-point maps from the
parser outputs.

RINEX 3 and RINEX 4 also support generating elevation and azimuth files. Differential
code bias (DCB) corrections from external files are not currently applied (WIP).

NB: The original code was written in 2015 with Python 2.7 and has since been updated to Python 3.9+ with a
focus on RINEX 3 / RINEX 4 support. The RINEX 2 parser is included for legacy support. The revamped new code has
been tested against the original teqc output with 2-decimal precision and is expected to produce equivalent
results. Lastly, the code has been refactored to use a more object-oriented approach, with a focus on maintainability and extensibility.
Copilot is used to generate the code structure and boilerplate, which has been modified to suit the specific needs of this project.


# RINEX 2 Parser

The current `rinex_2_tool.py` parses RINEX 2 observation files, prints a summary,
and writes the observations CSV to `Data_2/output/`.

Run it on the included default sample:

```sh
python3 rinex_2_tool.py
```

Specify a RINEX 2 observation file (place RINEX 2 files in the `Data/` directory):

```sh
python3 rinex_2_tool.py Data/kepa3150.13o
```

Recommended with `uv`, because you don’t need to activate the environment to run a command—use uv run, for example:

```shell
uv run python3 rinex_2_tool.py Data/kepa3150.13o
```

The parser API also returns a pandas DataFrame with one row per epoch/satellite.
Its index is `(epoch, satellite)` and its MultiIndex columns contain each
observation type's `value`, `lli` (Loss of Lock Indicator), and `signal_strength` fields.

The CSV is written automatically to `Data_2/output/`, using the input filename
with `.csv` appended. `--csv` changes the output filename while keeping it in
that folder:

```sh
uv run python3 rinex_2_tool.py Data_2/kepa3150.13o
uv run python3 rinex_2_tool.py Data_2/kepa3150.13o --csv Data_2/output/observations.csv
```

RINEX 2 navigation files can also be parsed; their CSV is written to the same
folder by default:

```sh
uv run python3 rinex_2_tool.py nav Data_2/brdc3170.13n
```

## Azimuth and Elevation
Generate compact azimuth and elevation files from an observation/navigation
pair (15-second sampling by default). These files are written to
`Data_2/output/` and use the compact format consumed by the legacy VTEC code:

```sh
uv run python3 rinex_2_tool.py angles Data_2/kepa3170.13o Data_2/brdc3170.13n
```

Use `--interval SECONDS` to change the sampling interval or `--outdir DIR` to
choose another output directory. Angle files cover the observation time span
available in the input file. The generated files were benchmarked against the
provided TEQC output over its overlapping two-hour span.

Calculate slant and vertical TEC from dual-frequency code and carrier-phase
observations plus the matching compact elevation file:

```sh
uv run python3 rinex_2_tool.py tec Data_2/kepa3170.13o
```

The CSV is written to `Data_2/output/kepa3170_tec.csv`. Code pseudoranges
(P1/P2 when available, otherwise C1/C2) provide the TEC level, while changes
in the L1/L2 carrier-phase geometry-free measurement smooth it with a
100-second Hatch filter. The output also includes a carrier-phase-only STEC
estimate, which is precise in relative changes but has an unresolved ambiguity
offset. A 350 km thin-shell mapping height is used by default
(`--shell-height KM` changes it). Since no DCB correction is applied, the
STEC/VTEC estimates are uncalibrated. Plotting is a separate step, so the TEC
CSV can be processed independently of the parser:

```sh
uv run python3 rinex_visualization.py tec Data_2/output/kepa3170_tec.csv
```

This writes per-satellite plots of phase-smoothed VTEC against a third-degree
polynomial fit, with residuals below, to
`Data_2/output/kepa3170_vtec_plots/`. `--ele FILE` selects another compact
elevation file and `--out FILE` selects the CSV destination. GPS-time epochs are
converted to UTC in the CSV and plot axis using the applicable leap-second
offset.

The reusable calculation API is in `tec_calculation.py`. It accepts normalized
observation rows (`epoch`, `satellite`, two code measurements, two phase
measurements, their frequencies, LLI flags, and `code_pair`) plus elevation
rows (`epoch`, `satellite`, `elevation_deg`). 

Signal selection and carrier frequencies are explicit so callers can choose
compatible signals for each constellation. Epochs in the observation and
elevation tables must use the same timescale.

Generate a pierce-point map from compact azimuth/elevation results using the
receiver's geodetic coordinates:

```sh
uv run python3 rinex_visualization.py map \
  --azimuth Data_2/output/kepa3170.azi \
  --elevation Data_2/output/kepa3170.ele \
  --latitude LATITUDE --longitude LONGITUDE
```

For a station-centered map independent of the RINEX parsers and satellite
inputs, use `map_generation.py` directly:

```sh
uv run python3 map_generation.py \
  --latitude 36.2 --longitude 34.3 --output output/station_map.png
```

An optional CSV with `elevation_deg` and `azimuth_deg` columns adds thin-shell
ionospheric pierce points:

```sh
uv run python3 map_generation.py \
  --latitude 36.2 --longitude 34.3 \
  --angles satellite_angles.csv \
  --shell-height 350 \
  --output output/station_ipp_map.png
```

The original `plot_map_globals`, `funky_earth_plot`, and `plot_map` plot types
are also available as independent functions and CLI modes. Global and globe
views need only station coordinates:

```sh
uv run python3 map_generation.py --mode global \
  --latitude LAT --longitude LON --output output/global.png
uv run python3 map_generation.py --mode funky-earth \
  --latitude LAT --longitude LON --output output/globe.png
```

For `--mode tracks`, pass a CSV with `satellite`, `latitude_deg`, and
`longitude_deg` columns:

```sh
uv run python3 map_generation.py --mode tracks \
  --latitude LAT --longitude LON --tracks pierce_points.csv \
  --output output/tracks.png
```

A test suite is included in `tests/` with synthetic RINEX 2 files generated by `tests/make_sample.py` (simulated GPS+GLONASS observables, 114 sats, cycle slips, LLI/SSI flags); `example_out/` holds pre-generated samples and full CLI output.            
Parser tests:

```sh
uv run python3 -m unittest discover -s tests -v
```


# RINEX 3/ RINEX 4 Parser and GNSS Data Extraction

Pure-Python parser for RINEX 3.x and RINEX 4.00 observation and navigation files
(GPS, Galileo, BeiDou first-class; GLONASS/QZSS/NavIC/SBAS parsed generically),
with timescale conversion (GPST/UTC/TAI/BDT/GST), tidy pandas
extraction and summary statistics. Plotting is provided separately by
`rinex_visualization.py`.

Requires: Python >= 3.9, numpy, pandas, matplotlib (pytest for the tests). These are however optional dependencies; the parser itself has no external dependencies.

## CLI

    python rinex_3_4_tool.py obs <file.o> [-o outdir] [--systems GEC]
    python rinex_3_4_tool.py nav <file.n> [-o outdir] [--systems GEC]
    python rinex_3_4_tool.py angles <file.o> <file.n> [-o outdir] [--interval 15] [--systems GEC]
    python rinex_3_4_tool.py tec <file.o> [--ele file.ele] [-o outdir]

`--systems` accepts one or more RINEX system letters (default: `GEC`):

| Letter | GNSS |
|---|---|
| `G` | GPS |
| `E` | Galileo |
| `C` | BeiDou |
| `R` | GLONASS |
| `J` | QZSS |
| `I` | NavIC (IRNSS) |
| `S` | SBAS |

For example, `--systems GREC` selects GPS, GLONASS, Galileo, and BeiDou.

RINEX 3 / RINEX 4 source files belong in `Data_3_4/input/`. By default, generated files
go to `Data_3_4/output/v304/` or `Data_3_4/output/v400/`, according to the
input version. Use `-o` to choose a different output root; the version-specific
subdirectory is still created. Every generated filename gets the local run date
in `DDMMYY` format, such as `observations_long_181126.csv`.

`obs` writes `observations_long_DDMMYY.csv` (one row per satellite/epoch/signal,
with time in the file scale plus GPST/UTC/TAI), `phase_pseudorange_DDMMYY.csv`,
`observations_wide_DDMMYY.csv` (one row per sat/epoch, one column per code),
`snr_stats_DDMMYY.csv`, `visibility_DDMMYY.csv`, `observation_counts_DDMMYY.csv`,
and `summary_DDMMYY.json`. All output filenames include the same run date.
Generate its SNR and satellite visibility plots separately from the long CSV:

```sh
uv run python3 rinex_visualization.py observations \
  Data_3_4/output/v400/observations_long_DDMMYY.csv
```
## Azimuth and Elevation
The shared visualization module also accepts compact azimuth/elevation files
from any parser that produces that format. `nav` writes
`ephemeris_DDMMYY.csv` (Keplerian + clock parameters, TOC in the file's timescale and GPST), plus `sto_DDMMYY.csv` / `ion_DDMMYY.csv` for RINEX 4
system-time offset and ionosphere records.

`angles` uses the observation receiver coordinates and broadcast ephemerides
to write compact `.azi` and `.ele` files under the selected output root's
version directory (GPS, Galileo, and BeiDou broadcast orbits are supported):

```sh
uv run python3 rinex_3_4_tool.py angles \
  Data_3_4/input/station.24o Data_3_4/input/brdc.24n
```

The files use the same `COMPACT2` TEQC format as the RINEX 2 geometry output and can
be supplied to the shared TEC calculation.
Use `--systems E` for a Galileo-only navigation file such as
`MERS00TUR_R_20230370000_01D_EN.rnx`. Geometry is computed only for satellites
present in both the observation epoch and the selected navigation ephemerides,
so selecting the available constellation also avoids unnecessary orbit work (this significantally reduces processing time).

Calculate uncalibrated STEC/VTEC from the RINEX 3/ RINEX 4 observations and generated
elevations, then plot the resulting CSV separately:

```sh
 uv run python3 rinex_3_4_tool.py angles Data_3_4/input/MERS00TUR_R_20230370000_01D_30S_MO.rnx \ 
  Data_3_4/input/MERS00TUR_R_20230370000_01D_EN.rnx --systems E
uv run python rinex_3_4_tool.py tec Data_3_4/input/MERS00TUR_R_20230370000_01D_30S_MO.rnx --system E
uv run python3 rinex_visualization.py tec Data_3_4/output/v303/MERS00TUR_R_20230370000_01D_30S_MO_tec.csv 
```

NB: The `systems` option is used to specify which satellite systems to process, and is corresponding to the letter preceding the 'N' in the navigation file name. 
For example, `MERS00TUR_R_20230370000_01D_EN.rnx` is a Galileo navigation file, so `--systems E` is used to select it.

By default, `tec` reads the matching `<observation-stem>.ele` file from the
versioned output directory and writes `<observation-stem>_tec.csv` there.
`--ele FILE` selects another elevation file. `--out DIR` selects an output
root, while `--out FILE.csv` writes to an explicit CSV path. GPS is selected by
default with `C1C/C2W` code and `L1C/L2W` phase signals at L1/L2 frequencies.
`--system E` and `--system C` select Galileo and BeiDou presets respectively;
use `--code1`, `--code2`, `--phase1`, `--phase2`, `--frequency1`, and
`--frequency2` to match the signals and frequencies actually present in your
file. The selected signal types must exist in the observation header. TEC
epochs are converted to UTC in the CSV; the calculation is uncalibrated unless
valid DCB corrections are applied separately. `.azi` is not needed for VTEC.

With `uv`, you don’t need to activate the environment to run a command—use `uv run`, for example:

```shell
uv run python rinex_3_4_tool.py obs Data_3_4/input/sample_v304.24o --systems E
```

## API

    from rinex_3_4_tool import read_obs, read_nav, to_gpst, convert_time, gps_week_tow

    header, df = read_obs("Data_3_4/input/SITE001.24o")  # df: tidy DataFrame
    header.approx_xyz                        # antenna XYZ from the header
    header.approx_llh                         # derived geodetic position
    header.obs_types                          # {"G": ["C1C", ...], ...}
    to_gpst(t, "BDT"); convert_time(t, "UTC", "TAI")
    nav = read_nav("Data_3_4/input/SITE001.24n"); nav.ephemeris

## Time scales

GPST–UTC = TAI–UTC − 19 s using the full leap-second table (1980–today);
GPST = GST, BDT = GPST − 14 s. UTC conversions are correct across leap
seconds (1972+). `gps_week_tow` returns continuous GPS week (no 1024 rollover)
and time-of-week.

## Tests

    uv run python3 -m pytest tests

Test fixtures are synthetic RINEX files generated by `tests/make_sample.py`
(simulated GPS+Galileo+BeiDou observables, 114 sats, cycle slips, LLI/SSI
flags); `example_out/` holds pre-generated samples and full CLI output.