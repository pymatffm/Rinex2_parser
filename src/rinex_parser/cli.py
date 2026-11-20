"""Unified command-line entry point for RINEX processing tools."""
from __future__ import annotations

import argparse
import sys
from typing import List, Optional


def main(argv: Optional[List[str]] = None) -> int:
    """Dispatch a command to its parser or visualization implementation."""
    arguments = list(sys.argv[1:] if argv is None else argv)
    parser = argparse.ArgumentParser(
        prog="rinex",
        description="Parse, process, and visualize RINEX GNSS data.",
        epilog=(
            "Commands: rinex2, rinex34, plot, map. Run "
            "'rinex <command> --help' for command-specific options."
        ),
    )
    parser.add_argument(
        "command",
        nargs="?",
        choices=("rinex2", "rinex34", "plot", "map"),
    )
    if not arguments or arguments[0] in ("-h", "--help"):
        parser.print_help()
        return 0

    command, *command_arguments = arguments
    handlers = {
        "rinex2": ("rinex_parser.rinex2.tool", "main"),
        "rinex34": ("rinex_parser.rinex34.tool", "main"),
        "plot": ("rinex_parser.visualization.plots", "main"),
        "map": ("rinex_parser.visualization.maps", "main"),
    }
    if command not in handlers:
        parser.error("invalid choice: {!r}".format(command))

    module_name, function_name = handlers[command]
    from importlib import import_module

    handler = getattr(import_module(module_name), function_name)
    return handler(command_arguments)
