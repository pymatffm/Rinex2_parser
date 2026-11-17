# Rinex2_parser

The current `main.py` parses RINEX 2 observation files and prints a summary; it
does not run the legacy navigation/VTEC processing pipeline.

Run it on the included sample:

```sh
python3 main.py
```

Or specify a RINEX 2 observation file:

```sh
python3 main.py Data/kepa3150.13o
```

The parser API also returns a pandas DataFrame with one row per epoch/satellite.
Its index is `(epoch, satellite)` and its MultiIndex columns contain each
observation type's `value`, `lli`, and `signal_strength` fields:

```python
from pathlib import Path
from main import read_observation_dataframe

header, observations = read_observation_dataframe(Path("Data/kepa3150.13o"))
print(observations.loc[(slice(None), "R21"), ("L1", "value")])
```

Export the table to CSV from the command line:

```sh
python3 main.py Data/kepa3150.13o --csv observations.csv
```

Parser tests:

```sh
python3 -m unittest discover -s tests -v
```
