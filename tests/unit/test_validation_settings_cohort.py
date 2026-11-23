"""Tests for the cohort comparison of read-end settings on constructed inputs.

Covers scripts/emseq_bedgraph_to_methylkit.py (strand from the reference and
autosomal totals), scripts/emseq_validation_settings_compare.R (tables from
constructed per-library totals and tile tables), and the read-end flags of
scripts/emseq_validation_settings_cohort.sh against workflows/validation.smk.
"""
import csv
import gzip
import io
import re
import shutil
import subprocess
from pathlib import Path

import pytest

import emseq_bedgraph_to_methylkit as b2m

REPO = Path(__file__).resolve().parents[2]


def write_fasta(path: Path, contigs: dict, width: int = 10) -> None:
    """FASTA with lines of `width` bases and its .fai index."""
    fai, offset = [], 0
    with open(path, "w") as fh:
        for name, seq in contigs.items():
            header = f">{name}\n"
            fh.write(header)
            offset += len(header)
            fai.append(f"{name}\t{len(seq)}\t{offset}\t{width}\t{width + 1}\n")
            for i in range(0, len(seq), width):
                line = seq[i:i + width] + "\n"
                fh.write(line)
                offset += len(line)
    Path(str(path) + ".fai").write_text("".join(fai))


@pytest.fixture
def ref(tmp_path):
    # contig 1: C at 0-based 2, G at 3, A at 12 (second line); chrX: C at 0
    fa = tmp_path / "ref.fa"
    write_fasta(fa, {"1": "AACGTTTTTTAAATT", "X": "CGAA"})
    return fa


def test_strand_from_reference_and_autosomal_totals(ref):
    bedgraph = io.StringIO(
        "track type=bedGraph\n"
        "1\t2\t3\t75\t3\t1\n"   # C: top strand
        "1\t3\t4\t50\t1\t1\n"   # G: bottom strand
        "1\t12\t13\t100\t2\t0\n"  # A: neither, counted as other_base
        "X\t0\t1\t100\t5\t0\n"  # sex chromosome: written, not in the totals
    )
    out = io.StringIO()
    totals = b2m.convert(bedgraph, str(ref), out)
    rows = out.getvalue().splitlines()
    assert rows[0] == "chrBase\tchr\tbase\tstrand\tcoverage\tfreqC\tfreqT"
    assert rows[1] == "1.3\t1\t3\tF\t4\t75.00\t25.00"
    assert rows[2] == "1.4\t1\t4\tR\t2\t50.00\t50.00"
    assert rows[3] == "X.1\tX\t1\tF\t5\t100.00\t0.00"
    assert totals == {"positions": 2, "meth": 4, "unmeth": 2, "other_base": 1}


def test_contig_sequence_spans_line_breaks(ref):
    index = b2m.read_fai(str(ref) + ".fai")
    with open(ref, "rb") as fh:
        assert b2m.contig_sequence(fh, index["1"]) == b"AACGTTTTTTAAATT"


def tile_table(path: Path, rows: list) -> None:
    with gzip.open(path, "wt") as fh:
        fh.write("#Date:2026-11-22\n#Class:methylDiffDB\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript not on PATH")
def test_compare_tables(tmp_path):
    if subprocess.run(["Rscript", "-e", "library(data.table)"], capture_output=True).returncode:
        pytest.skip("data.table not installed")
    work = tmp_path / "work"
    (work / "global").mkdir(parents=True)
    # four libraries; neb adds 1 point to healthy and 0.5 to patients, all10 adds 2 and 1
    base = {"NH_1": 80, "NH_2": 78, "PRO_1": 74, "PRO_2": 72}
    add = {"none": (0, 0), "neb": (1, 0.5), "all10": (2, 1)}
    for lib, pct in base.items():
        for s, (healthy, patient) in add.items():
            # the second library of each group gains 0.1 point more under a mask
            spread = 0.1 if lib.endswith("2") and s != "none" else 0
            p = pct + (healthy if lib.startswith("NH_") else patient) + spread
            meth = int(round(p * 100))
            (work / "global" / f"{lib}.{s}.txt").write_text(
                f"positions=10\tmeth={meth}\tunmeth={10000 - meth}\tother_base=0\n")
    tiles = {
        "none": [[1, 1, 10000, "*", 1e-5, 1e-3, -30], [1, 10001, 20000, "*", 0.5, 0.6, 5]],
        "neb": [[1, 1, 10000, "*", 1e-5, 1e-3, -31], [1, 10001, 20000, "*", 1e-4, 1e-3, 26]],
        "all10": [[1, 1, 10000, "*", 1e-5, 1e-3, 30], [1, 10001, 20000, "*", 1e-4, 1e-3, 27]],
    }
    for s, rows in tiles.items():
        (work / "dmr" / s).mkdir(parents=True)
        tile_table(work / "dmr" / s / "methylDiff_x.tiled.txt.bgz", rows)
    index_map = tmp_path / "index.tsv"
    index_map.write_text("# shuffled\nlibrary\tlibrary_index\nPRO_2\t1\nNH_1\t2\nPRO_1\t3\nNH_2\t4\n")
    out = tmp_path / "out"
    subprocess.run(["Rscript", str(REPO / "scripts" / "emseq_validation_settings_compare.R"), str(work), str(out),
                    str(index_map)], check=True, capture_output=True)

    glob = list(csv.DictReader(open(out / "global.tsv"), delimiter="\t"))
    assert "library" not in glob[0]
    assert [r["library_index"] for r in glob[:3]] == ["1", "1", "1"]
    assert {r["group"] for r in glob if r["library_index"] == "2"} == {"healthy"}

    diff = {r["comparison"]: r for r in csv.DictReader(open(out / "global_diff.tsv"), delimiter="\t")}
    assert float(diff["all10 minus none"]["estimate"]) == pytest.approx(1.55)
    assert float(diff["patient minus healthy, of neb minus none"]["estimate"]) == pytest.approx(-0.5)
    assert float(diff["patient minus healthy, none"]["estimate"]) == pytest.approx(-6.0)

    called = {r["setting"]: r for r in csv.DictReader(open(out / "tiles.tsv"), delimiter="\t")}
    assert (called["none"]["called_q01_d25"], called["neb"]["called_q01_d25"]) == ("1", "2")

    pairs = list(csv.DictReader(open(out / "tiles_pair.tsv"), delimiter="\t"))
    strict = {r["comparison"]: r for r in pairs if r["q"] == "0.01"}
    assert strict["neb against none"]["called_both"] == "1"
    assert float(strict["neb against none"]["jaccard"]) == pytest.approx(0.5)
    # the first tile is called under none and all10 with opposite signs
    assert strict["all10 against none"]["sign_differs"] == "1"


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript not on PATH")
def test_compare_stops_on_library_missing_from_map(tmp_path):
    if subprocess.run(["Rscript", "-e", "library(data.table)"], capture_output=True).returncode:
        pytest.skip("data.table not installed")
    (tmp_path / "global").mkdir()
    (tmp_path / "global" / "NH_9.none.txt").write_text("positions=1\tmeth=1\tunmeth=1\tother_base=0\n")
    index_map = tmp_path / "index.tsv"
    index_map.write_text("library\tlibrary_index\nNH_1\t1\n")
    run = subprocess.run(["Rscript", str(REPO / "scripts" / "emseq_validation_settings_compare.R"), str(tmp_path),
                          str(tmp_path / "out"), str(index_map)], capture_output=True, text=True)
    assert run.returncode != 0
    assert "no entry for: NH_9" in run.stderr


def test_driver_flags_match_workflow():
    smk = (REPO / "workflows" / "validation.smk").read_text()
    match = re.search(r"VALIDATION_END_SETTINGS = \{(.*?)\}", smk, re.S)
    assert match
    block = match.group(1)
    workflow = dict(re.findall(r'"(\w+)": "([^"]*)"', block))
    driver = (REPO / "scripts" / "emseq_validation_settings_cohort.sh").read_text()
    match = re.search(r"declare -A FLAGS=\((.*)\)", driver)
    assert match
    line = match.group(1)
    flags = dict(re.findall(r'\[(\w+)\]="([^"]*)"', line))
    assert flags == workflow
