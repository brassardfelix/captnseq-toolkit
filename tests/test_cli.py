import subprocess
import sys

import pytest

from captnseq_toolkit.cli import SUBCOMMANDS, main


def test_help_lists_all_subcommands(capsys):
    assert main(["--help"]) == 0
    out = capsys.readouterr().out
    for name in SUBCOMMANDS:
        assert name in out


def test_unknown_subcommand():
    assert main(["nope"]) == 2


@pytest.mark.parametrize("name", list(SUBCOMMANDS))
def test_subcommand_help(name):
    r = subprocess.run([sys.executable, "-m", "captnseq_toolkit", name, "--help"], capture_output=True, text=True)
    assert r.returncode == 0
    assert f"captnseq-toolkit {name}" in r.stdout
