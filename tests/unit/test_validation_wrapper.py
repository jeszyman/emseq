# tests/unit/test_validation_wrapper.py
"""Structure of the validation module and its wrapper.

The module (workflows/emseq_validation.smk) reads no config, and every
MethylDackel command in it keeps pairs that lack the proper-pair flag. The
wrapper (workflows/validation.smk) takes its config on the command line and
stops when a library has no index for the result files.
"""
import re
import shutil
import subprocess

import pytest
import yaml

from conftest import REPO_ROOT

MODULE = REPO_ROOT / "workflows" / "emseq_validation.smk"
WRAPPER = REPO_ROOT / "workflows" / "validation.smk"


def methyldackel_commands(text):
    """Shell commands that start a line with MethylDackel, with their continuation lines."""
    return [m.strip() for m in re.findall(r"^[ \t]+MethylDackel (?:extract|mbias)(?:[^\n]*\\\n)*[^\n]*", text, re.M)]


def test_module_reads_no_config():
    text = MODULE.read_text()
    assert "config[" not in text
    assert "configfile:" not in text
    assert not re.search(r"^D_\w+\s*=", text, re.M)


def test_wrapper_names_no_config_file():
    assert "configfile:" not in WRAPPER.read_text()


def test_every_methyldackel_command_keeps_unflagged_pairs():
    commands = methyldackel_commands(MODULE.read_text())
    assert commands
    for command in commands:
        assert "--keepDiscordant" in command


def test_mbias_prints_suggested_bounds_and_lifts_the_depth_cap():
    mbias = [c for c in methyldackel_commands(MODULE.read_text()) if c.startswith("MethylDackel mbias")]
    assert mbias
    for command in mbias:
        assert "--noSVG" not in command
        assert re.search(r"-D \d+", command)


def base_config(tmp_path, index):
    return {
        "main-data-dir": str(tmp_path / "data"),
        "envs": {"emseq": "config/emseq-conda-env.yaml", "methylkit": "config/methylkit-conda-env.yaml"},
        "repos": {"emseq": str(REPO_ROOT)},
        "library-ids": ["libA", "libB"],
        "keep-bed": "keep.bed",
        "exclude-bed": "exclude.bed.gz",
        "emseq_ref_assemblies": {
            "host": {"name": "host", "input": "host.fa.gz"},
            "puc19": {"name": "puc19", "input": "pUC19.fa.gz"},
            "unmeth_lambda": {"name": "unmeth_lambda", "input": "lambda.fa.gz"},
        },
        "validation-host-ref": "host",
        "validation-dataset": "unit-test",
        "validation-library-index": index,
    }


@pytest.mark.skipif(shutil.which("snakemake") is None, reason="snakemake not on PATH")
def test_wrapper_stops_and_names_a_library_without_an_index(tmp_path):
    config = tmp_path / "config.yaml"
    config.write_text(yaml.safe_dump(base_config(tmp_path, {"libA": 1})))
    proc = subprocess.run(
        ["snakemake", "-s", str(WRAPPER), "--configfile", str(config), "--cores", "1", "-n"],
        capture_output=True, text=True, cwd=REPO_ROOT,
    )
    assert proc.returncode != 0
    assert "validation-library-index has no entry for: libB" in proc.stdout + proc.stderr
