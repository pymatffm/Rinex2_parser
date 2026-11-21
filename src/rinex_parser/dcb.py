"""Parsing and application of CODE differential code-bias files."""
from __future__ import annotations

import re
from pathlib import Path
from typing import Dict

_SPEED_OF_LIGHT = 299792458.0
_SATELLITE_ROW = re.compile(
    r"^(?P<system>[GRECJIS])(?P<number>\d{2})\s+"
    r"(?P<value>[+-]?(?:\d+(?:\.\d*)?|\.\d+))\s+"
    r"(?P<rms>[+-]?(?:\d+(?:\.\d*)?|\.\d+))"
)


def read_dcb_file(path) -> Dict[str, float]:
    """Read satellite DCB values from a CODE monthly ASCII file.

    The returned mapping uses RINEX satellite IDs (for example ``G01``) and
    metre values. Station receiver rows in the CODE format are ignored.
    """
    path = Path(path)
    kind = ""
    satellites: Dict[str, float] = {}

    with path.open("r", encoding="ascii", errors="replace") as stream:
        for line in stream:
            upper = line.upper()
            if "P1-P2" in upper:
                kind = "P1-P2"
            elif "P1-C1" in upper:
                kind = "P1-C1"

            match = _SATELLITE_ROW.match(line)
            if not match:
                continue
            satellite = "{}{}".format(match.group("system"), match.group("number"))
            satellites[satellite] = (
                float(match.group("value")) * 1e-9 * _SPEED_OF_LIGHT
            )

    if not kind:
        raise ValueError("{} does not identify a supported CODE DCB type".format(path))
    if not satellites:
        raise ValueError("{} contains no satellite DCB values".format(path))
    return satellites


def geometry_free_code_correction(
    satellite: str,
    code_1: str,
    code_2: str,
    p1p2_m: Dict[str, float],
    p1c1_m: Dict[str, float],
) -> float:
    """Return metres to add to ``code_2 - code_1``.

    P1-P2 and P1-C1 corrections are kept separate because C1 observables are
    not necessarily P1 observables. The function supports C1/P1 and C2/P2
    geometry-free pairs; unsupported combinations are left uncorrected.
    """
    correction = 0.0
    if code_1 in {"C1C", "C1W", "C1X", "C1", "P1"}:
        if code_2 in {"C2W", "C2X", "C2", "P2"}:
            correction += p1p2_m.get(satellite, 0.0)
            if code_1.startswith("C") and code_1 != "P1":
                correction -= p1c1_m.get(satellite, 0.0)
    return correction
