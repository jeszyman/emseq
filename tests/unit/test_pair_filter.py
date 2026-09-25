# tests/unit/test_pair_filter.py
"""The BAM pair filter keeps dinucleosome-length fragments.

bwa-meth (BWA-MEM) sets the proper-pair flag (0x2) from the insert-size
distribution it infers per library, so on cfDNA it stops flagging pairs at
about 270-280 bp and a `-f 2` filter drops every dinucleosome-length fragment.
EMSEQ_PAIR_FILTER in workflows/emseq.smk replaces that flag with an explicit
check: same chromosome, inward-facing mates, |TLEN| <= 1000.
"""
import re
import shutil
import subprocess

import pytest

from conftest import REPO_ROOT

SMK = REPO_ROOT / "workflows" / "emseq.smk"

HEADER = "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:100000\n@SQ\tSN:chr2\tLN:100000\n"
SEQ = "A" * 50
QUAL = "I" * 50


def pair(name, flag1, flag2, pos1, pos2, tlen, rnext="="):
    """Two SAM records for one pair; flags exclude 0x2 unless given."""
    r1 = f"{name}\t{flag1}\tchr1\t{pos1}\t60\t50M\t{rnext}\t{pos2}\t{tlen}\t{SEQ}\t{QUAL}\n"
    r2 = f"{name}\t{flag2}\tchr1\t{pos2}\t60\t50M\t{rnext}\t{pos1}\t{-tlen}\t{SEQ}\t{QUAL}\n"
    return r1 + r2


# 0x1 paired, 0x40/0x80 read 1/2, 0x10 reverse, 0x20 mate reverse
RECORDS = {
    # inward-facing, 150 bp, flagged proper by BWA
    "mono_proper": pair("mono_proper", 1 | 2 | 0x40 | 0x20, 1 | 2 | 0x80 | 0x10, 1000, 1101, 150),
    # inward-facing, 330 bp, NOT flagged proper (the case -f 2 drops)
    "di_unflagged": pair("di_unflagged", 1 | 0x40 | 0x20, 1 | 0x80 | 0x10, 2000, 2281, 330),
    # inward-facing, 900 bp, not flagged proper
    "long_unflagged": pair("long_unflagged", 1 | 0x40 | 0x20, 1 | 0x80 | 0x10, 5000, 5851, 900),
    # inward-facing but 5000 bp: likely mismapped or structural
    "too_long": pair("too_long", 1 | 0x40 | 0x20, 1 | 0x80 | 0x10, 10000, 14951, 5000),
    # both mates on the same strand
    "same_strand": pair("same_strand", 1 | 0x40, 1 | 0x80, 20000, 20201, 250),
    # outward-facing: reverse mate leftmost (positive TLEN, per the SAM spec), forward mate downstream
    "outward": pair("outward", 1 | 0x40 | 0x10, 1 | 0x80 | 0x20, 30000, 30201, 250),
}
MATE_OTHER_CHROM = (
    f"trans\t{1 | 0x40 | 0x20}\tchr1\t40000\t60\t50M\tchr2\t40000\t0\t{SEQ}\t{QUAL}\n"
    f"trans\t{1 | 0x80 | 0x10}\tchr2\t40000\t60\t50M\tchr1\t40000\t0\t{SEQ}\t{QUAL}\n"
)
KEEP = {"mono_proper", "di_unflagged", "long_unflagged"}
DROP = {"too_long", "same_strand", "outward", "trans"}


def pair_filter_expr():
    m = re.search(r'^EMSEQ_PAIR_FILTER\s*=\s*\(?\s*((?:"[^"]*"\s*)+)\)?', SMK.read_text(), re.M)
    assert m, "EMSEQ_PAIR_FILTER is not defined in workflows/emseq.smk"
    return "".join(re.findall(r'"([^"]*)"', m.group(1)))


@pytest.mark.skipif(shutil.which("samtools") is None, reason="samtools not on PATH")
def test_pair_filter_keeps_long_fragments_and_drops_bad_pairs(tmp_path):
    sam = tmp_path / "pairs.sam"
    sam.write_text(HEADER + "".join(RECORDS.values()) + MATE_OTHER_CHROM)
    out = subprocess.run(
        ["samtools", "view", "-f", "1", "-F", "12", "-e", pair_filter_expr(), str(sam)],
        check=True, capture_output=True, text=True,
    ).stdout
    kept = [line.split("\t")[0] for line in out.splitlines()]
    assert set(kept) == KEEP
    # both mates of every kept pair survive, so dupsifter never sees an orphan
    assert all(kept.count(name) == 2 for name in KEEP)


def rule_shell(rule):
    text = SMK.read_text()
    start = text.index(f"rule {rule}:")
    end = text.find("\nrule ", start + 1)
    return text[start:end if end != -1 else None]


@pytest.mark.parametrize("rule", ["emseq_dedup", "emseq_filter_bam"])
def test_rules_use_pair_filter_not_proper_pair_flag(rule):
    body = rule_shell(rule)
    assert not re.search(r"samtools view[^\n]*-f\s*(0x)?2\b", body), f"{rule} still filters on the proper-pair flag"
    assert "EMSEQ_PAIR_FILTER" in body, f"{rule} does not apply EMSEQ_PAIR_FILTER"
