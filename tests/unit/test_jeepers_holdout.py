"""Opposite-strand pairs and the hold-out BAMs of emseq_jeepers_holdout.py."""
import pathlib
import sys

import pytest

pysam = pytest.importorskip("pysam")

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import emseq_jeepers_holdout as ho  # noqa: E402
from jeepers_fixtures import reads_of, sam_pair, write_bam  # noqa: E402


def test_same_coordinates_opposite_strands_pair(tmp_path):
    bam = write_bam(tmp_path / "a.bam", [sam_pair("t", 101, 300), sam_pair("b", 101, 300, strand="OB")])
    assert ho.strand_pairs(bam) == [("t", "b")]
    assert ho.near_pairs(bam) == []


def test_same_strand_same_coordinates_do_not_pair(tmp_path):
    bam = write_bam(tmp_path / "a.bam", [sam_pair("t1", 101, 300), sam_pair("t2", 101, 300)])
    assert ho.strand_pairs(bam) == []
    # with several fragments on a strand at the same coordinates, the first is taken
    bam = write_bam(
        tmp_path / "b.bam",
        [sam_pair("t1", 101, 300), sam_pair("t2", 101, 300), sam_pair("b", 101, 300, strand="OB")],
    )
    assert ho.strand_pairs(bam) == [("t1", "b")]


def test_end_off_by_one_is_near_not_same(tmp_path):
    bam = write_bam(
        tmp_path / "a.bam",
        [
            sam_pair("t", 101, 300),
            sam_pair("end_plus", 101, 301, strand="OB"),
            sam_pair("start_minus", 100, 300, strand="OB"),
            sam_pair("both_off", 102, 301, strand="OB"),
            sam_pair("off_by_two", 101, 302, strand="OB"),
        ],
    )
    assert ho.strand_pairs(bam) == []
    assert sorted(ho.near_pairs(bam)) == [("t", "end_plus"), ("t", "start_minus")]


PAIRS = [(f"t{i}", f"b{i}") for i in range(200)]


def test_split_is_stable_for_the_seed():
    assert ho.split(PAIRS) == ho.split(PAIRS) == ho.split(PAIRS, seed=20261120)
    assert ho.split(PAIRS) != ho.split(PAIRS, seed=1)


def test_split_puts_each_pair_in_both_sets_once():
    a, b = ho.split(PAIRS)
    assert not a & b
    for top, bottom in PAIRS:
        assert (top in a and bottom in b) or (top in b and bottom in a)
    tops_in_a = sum(top in a for top, _ in PAIRS)
    assert 60 < tops_in_a < 140


def test_out_a_lacks_both_reads_of_set_a_fragments_and_out_b_has_them(tmp_path):
    pairs = []
    for i in range(20):
        start = 101 + 200 * i
        pairs += [sam_pair(f"t{i}", start, start + 180), sam_pair(f"b{i}", start, start + 180, strand="OB")]
    pairs.append(sam_pair("alone", 4500, 4700))
    bam = write_bam(tmp_path / "in.bam", pairs)
    out_a, out_b, table = tmp_path / "a.bam", tmp_path / "b.bam", tmp_path / "pairs.tsv"
    ho.main(["--bam", bam, "--out-a", str(out_a), "--out-b", str(out_b), "--out-pairs", str(table)])
    set_a, set_b = ho.split(ho.strand_pairs(bam))
    assert len(set_a) == len(set_b) == 20
    names_a = [r.query_name for r in reads_of(out_a)]
    names_b = [r.query_name for r in reads_of(out_b)]
    assert not set(names_a) & set_a
    assert not set(names_b) & set_b
    for name in set_a:
        assert names_b.count(name) == 2
    for name in set_b:
        assert names_a.count(name) == 2
    assert names_a.count("alone") == names_b.count("alone") == 2
    rows = [line.split("\t") for line in table.read_text().rstrip("\n").split("\n")[1:]]
    assert len(rows) == 20
    for kind, top, bottom, removed_in in rows:
        assert kind == "same"
        assert removed_in == ("A" if top in set_a else "B")
        assert (bottom in set_b) == (top in set_a)
