import subprocess
import sys

import pandas as pd
import pytest

from captnseq_toolkit import __version__
from captnseq_toolkit.cli import SUBCOMMANDS, main


def test_help_lists_all_subcommands(capsys):
    with pytest.raises(SystemExit) as exc:
        main(["--help"])
    assert exc.value.code == 0
    out = capsys.readouterr().out
    for name in SUBCOMMANDS:
        assert name in out


def test_version(capsys):
    with pytest.raises(SystemExit) as exc:
        main(["--version"])
    assert exc.value.code == 0
    assert capsys.readouterr().out.strip() == f"captnseq-toolkit {__version__}"


def test_no_subcommand_is_an_error(capsys):
    with pytest.raises(SystemExit) as exc:
        main([])
    assert exc.value.code == 2
    assert "required" in capsys.readouterr().err


def test_unknown_subcommand_is_an_error(capsys):
    with pytest.raises(SystemExit) as exc:
        main(["nope"])
    assert exc.value.code == 2
    assert "invalid choice" in capsys.readouterr().err


@pytest.mark.parametrize("name", list(SUBCOMMANDS))
def test_subcommand_help(name):
    r = subprocess.run([sys.executable, "-m", "captnseq_toolkit", name, "--help"], capture_output=True, text=True)
    assert r.returncode == 0
    assert f"captnseq-toolkit {name}" in r.stdout


@pytest.mark.parametrize("name", list(SUBCOMMANDS))
def test_subcommand_missing_required_options(name, capsys):
    with pytest.raises(SystemExit) as exc:
        main([name])
    assert exc.value.code == 2
    assert "required" in capsys.readouterr().err


def test_sites2genes_runs_through_top_level_cli(tmp_path):
    features = tmp_path / "cds.bed"
    features.write_text("chr\t0\t1000\tgeneA\t0\t1\n")
    sample = tmp_path / "intersect.bed"
    sample.write_text("chr\t0\t1000\tgeneA\t0\t1\tchr\t500\t501\t10\n")
    out = tmp_path / "out.tsv"
    rc = main([
        "sites2genes",
        "--features", str(features),
        "--sample", str(sample),
        "--length-pc-start", "0",
        "--length-pc-stop", "100",
        "--output", str(out),
    ])
    assert rc == 0
    df = pd.read_csv(out, sep="\t").set_index("gene_name")
    assert df.loc["geneA", "ins_count"] == 1
