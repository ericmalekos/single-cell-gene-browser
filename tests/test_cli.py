"""The scbrowser command line."""
import pytest

from scbrowser.cli import main


def test_check_runs(capsys):
    assert main(["check"]) == 0
    out = capsys.readouterr().out
    assert "Inputs" in out and "Stores" in out


def test_version(capsys):
    with pytest.raises(SystemExit) as exc:
        main(["--version"])
    assert exc.value.code == 0 and "scbrowser" in capsys.readouterr().out


def test_unknown_options_only_for_serve():
    with pytest.raises(SystemExit):
        main(["check", "--bogus"])
