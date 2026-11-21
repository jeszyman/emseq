# tests/unit/test_validation_fragment_end.py
"""Known-answer tests for scripts/emseq_validation_fragment_end.py.

The script tallies CpG calls by distance from the two ends of the original
strand. End repair fills the 3' end of each strand with unmethylated
cytosine, so the 3' tally shows the fill-in loss and the 5' tally is its
control. Pairs are built by hand on a reference that is all A except for
CpGs at chosen positions.

Coordinates in comments are 0-based, end-exclusive.
"""
import shutil
import subprocess

import pytest

pysam = pytest.importorskip("pysam")

import emseq_validation_fragment_end as fe

pytestmark = pytest.mark.skipif(shutil.which("samtools") is None, reason="samtools not on PATH")

# 0x1 paired, 0x40/0x80 read 1/2, 0x10 reverse, 0x20 mate reverse; no 0x2
OT_R1, OT_R2 = 1 | 0x40 | 0x20, 1 | 0x80 | 0x10
OB_R1, OB_R2 = 1 | 0x40 | 0x10, 1 | 0x80 | 0x20

# CpG positions (C at the first number)
CPGS = [50, 197, 251, 897]


def reference():
    seq = ["A"] * 1000
    for c in CPGS:
        seq[c], seq[c + 1] = "C", "G"
    return "".join(seq)


REF = reference()


def read_seq(start, length, changes=None):
    seq = list(REF[start:start + length])
    for pos, base in (changes or {}).items():
        seq[pos - start] = base
    return "".join(seq)


def record(name, flag, start, length, mate_start, tlen, changes=None):
    seq = read_seq(start, length, changes)
    return (f"{name}\t{flag}\tctg\t{start + 1}\t60\t{length}M\t=\t{mate_start + 1}\t{tlen}\t"
            f"{seq}\t{'I' * length}\tMC:Z:{length}M\n")


PAIRS = {
    # OT fragment 100-200 read through by both mates; methylated CpG at 197,
    # 3 bases from the strand's 3' end and 98 from its 5' end
    "ot_full": record("ot_full", OT_R1, 100, 100, 100, 100) + record("ot_full", OT_R2, 100, 100, 100, -100),
    # OB fragment 250-350 read through by both mates; unmethylated CpG whose G
    # is at 252, 3 bases from the strand's 3' end (the left end of the fragment)
    "ob_full": record("ob_full", OB_R2, 250, 100, 250, 100, {252: "A"})
    + record("ob_full", OB_R1, 250, 100, 250, -100, {252: "A"}),
    # OT fragment 0-100; the mates disagree at the CpG at 50
    "ot_disagree": record("ot_disagree", OT_R1, 0, 100, 0, 100)
    + record("ot_disagree", OT_R2, 0, 100, 0, -100, {50: "T"}),
    # OT fragment 600-900, mates apart; methylated CpG at 897 is covered by
    # read 2 alone
    "ot_apart": record("ot_apart", OT_R1, 600, 50, 850, 300) + record("ot_apart", OT_R2, 850, 50, 600, -300),
}
HEADER = "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:ctg\tLN:1000\n"


@pytest.fixture(scope="module")
def workdir(tmp_path_factory):
    d = tmp_path_factory.mktemp("fragment_end")
    (d / "ref.fa").write_text(f">ctg\n{REF}\n")
    subprocess.run(["samtools", "faidx", str(d / "ref.fa")], check=True)
    return d


def tally(workdir, names, **kwargs):
    sam = workdir / "pairs.sam"
    sam.write_text(HEADER + "".join(PAIRS[n] for n in names))
    bam = workdir / "pairs.bam"
    subprocess.run(["samtools", "sort", "-o", str(bam), str(sam)], check=True, capture_output=True)
    subprocess.run(["samtools", "index", str(bam)], check=True)
    return fe.tally(bam, workdir / "ref.fa", **kwargs)


def calls(ends, key):
    """(n_methylated, n_unmethylated) for one (end, mate, distance) key."""
    row = ends.get(key)
    return None if row is None else (row["n_methylated"], row["n_unmethylated"])


def test_read_through_pair_reports_the_three_prime_cpg_from_both_mates(workdir):
    ends, _ = tally(workdir, ["ot_full"])
    assert calls(ends, ("3p", 2, 3)) == (1, 0)
    assert calls(ends, ("3p", 1, 3)) == (1, 0)
    assert calls(ends, ("5p", 1, 98)) == (1, 0)
    assert calls(ends, ("5p", 2, 98)) == (1, 0)


def test_overlapped_cpg_is_counted_once_per_fragment_in_the_pooled_rows(workdir):
    ends, _ = tally(workdir, ["ot_full"])
    assert calls(ends, ("3p", 0, 3)) == (1, 0)
    assert calls(ends, ("5p", 0, 98)) == (1, 0)


def test_bottom_strand_three_prime_end_is_the_left_end_of_the_fragment(workdir):
    ends, _ = tally(workdir, ["ob_full"])
    assert calls(ends, ("3p", 2, 3)) == (0, 1)
    assert calls(ends, ("3p", 0, 3)) == (0, 1)
    assert calls(ends, ("5p", 0, 98)) == (0, 1)


def test_mates_apart_report_the_three_prime_cpg_from_read_2_alone(workdir):
    ends, _ = tally(workdir, ["ot_apart"])
    assert calls(ends, ("3p", 2, 3)) == (1, 0)
    assert calls(ends, ("3p", 1, 3)) is None
    assert calls(ends, ("3p", 0, 3)) == (1, 0)


def test_pair_without_the_proper_pair_flag_is_counted(workdir):
    # none of the constructed pairs carries flag 0x2
    ends, _ = tally(workdir, ["ot_apart"])
    assert sum(r["n_methylated"] + r["n_unmethylated"] for (e, m, d), r in ends.items() if m == 0 and e == "3p") == 1


def test_distance_beyond_the_maximum_is_left_out(workdir):
    ends, _ = tally(workdir, ["ot_full"], max_distance=10)
    assert calls(ends, ("3p", 0, 3)) == (1, 0)
    assert calls(ends, ("5p", 0, 98)) is None


def test_disagreeing_mates_are_tallied_at_both_read_positions(workdir):
    ends, agreement = tally(workdir, ["ot_disagree"])
    # CpG at 50 of fragment 0-100: read 1 cycle 51, read 2 (reverse) cycle 50
    assert agreement[(1, 51)] == {"n_agree": 0, "n_disagree": 1}
    assert agreement[(2, 50)] == {"n_agree": 0, "n_disagree": 1}
    # the pooled row takes read 1's call
    assert calls(ends, ("5p", 0, 51)) == (1, 0)


def test_agreeing_mates_are_tallied(workdir):
    _, agreement = tally(workdir, ["ot_full"])
    # CpG at 197 of fragment 100-200: read 1 cycle 98, read 2 (reverse) cycle 3
    assert agreement[(1, 98)] == {"n_agree": 1, "n_disagree": 0}
    assert agreement[(2, 3)] == {"n_agree": 1, "n_disagree": 0}


def test_neb_mask_removes_the_call_only_where_read_1_does_not_cover_it(workdir):
    # read 2 cycle 3 is inside the NEB mask (4 bases on OT pairs)
    apart, _ = tally(workdir, ["ot_apart"])
    assert apart[("3p", 0, 3)]["n_kept_neb"] == 0
    full, _ = tally(workdir, ["ot_full"])
    assert full[("3p", 0, 3)]["n_kept_neb"] == 1


def test_ten_base_mask_removes_the_call_from_both_mates(workdir):
    # read 1 cycle 98 of 100 and read 2 cycle 3 are both within 10 of a read end
    full, _ = tally(workdir, ["ot_full"])
    assert full[("3p", 0, 3)]["n_kept_all10"] == 0
    # the CpG at 50 of fragment 0-100 is mid-read in both mates
    mid, _ = tally(workdir, ["ot_disagree"])
    assert mid[("5p", 0, 51)]["n_kept_all10"] == 1


def test_minimum_fragment_length_leaves_out_shorter_fragments(workdir):
    # ot_full is a 100 bp fragment and ot_apart a 300 bp fragment
    ends, _ = tally(workdir, ["ot_full", "ot_apart"], min_fragment=160)
    assert calls(ends, ("3p", 0, 3)) == (1, 0)
    assert calls(ends, ("3p", 1, 3)) is None


def test_script_writes_both_fragment_sets(workdir, tmp_path):
    import sys
    from conftest import REPO_ROOT
    sam = workdir / "pairs.sam"
    sam.write_text(HEADER + PAIRS["ot_full"] + PAIRS["ot_apart"])
    bam = workdir / "cli.bam"
    subprocess.run(["samtools", "sort", "-o", str(bam), str(sam)], check=True, capture_output=True)
    subprocess.run(["samtools", "index", str(bam)], check=True)
    out = tmp_path / "ends.tsv"
    subprocess.run([sys.executable, str(REPO_ROOT / "scripts" / "emseq_validation_fragment_end.py"),
                    "--bam", str(bam), "--fasta", str(workdir / "ref.fa"),
                    "--out-ends", str(out), "--out-agreement", str(tmp_path / "agree.tsv")], check=True)
    rows = [line.split("\t") for line in out.read_text().splitlines()]
    assert rows[0][:4] == ["fragments", "end", "mate", "distance"]
    pooled = {r[0]: int(r[4]) for r in rows[1:] if r[1:4] == ["3p", "0", "3"]}
    assert pooled == {"all": 2, "long": 1}
