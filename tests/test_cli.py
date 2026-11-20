import pytest

from rinex_parser.cli import main


@pytest.mark.parametrize(
    "command",
    ("rinex2", "rinex34", "plot", "map"),
)
def test_cli_dispatches_to_command_help(command, capsys):
    with pytest.raises(SystemExit) as exc_info:
        main([command, "--help"])

    assert exc_info.value.code == 0
    assert "usage:" in capsys.readouterr().out
