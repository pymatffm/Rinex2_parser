"""Compatibility wrapper for the reusable TEC calculation API."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))
from rinex_parser import tec_calculation as _implementation

sys.modules[__name__] = _implementation
