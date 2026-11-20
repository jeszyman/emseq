# tests/unit/test_methyldackel_unflagged_pairs.py
"""MethylDackel counts read pairs that lack the proper-pair flag.

MethylDackel ignores a paired alignment whose proper-pair flag (0x2) is unset
unless it is given --keepDiscordant. bwa-meth (BWA-MEM) sets that flag from an
insert-size range it infers per batch of reads, so it leaves the flag unset on
cfDNA fragments of about 280 bp or longer, and on nearly every pUC19 control
pair (too few pairs per batch to infer a range). The pipeline keeps those pairs
at filtering (EMSEQ_PAIR_FILTER), so each MethylDackel command must keep them
at calling too.
"""
import re
import shutil
import subprocess

import pytest

from conftest import REPO_ROOT

SMK = REPO_ROOT / "workflows" / "emseq.smk"

RULES = ["emseq_methyldackel_spike", "emseq_mbias", "emseq_methyldackel"]

# 300 bp reference with a CpG every 10 bases
REF = "ATCGATTAGC" * 30
QUAL = "I" * 50
# One inward-facing pair, 250 bp, MAPQ 60, proper-pair flag unset.
# 0x1 paired, 0x40/0x80 read 1/2, 0x10 reverse, 0x20 mate reverse.
# SEQ equals the reference, so every cytosine reads as methylated.
PAIR = (
    f"unflagged\t{1 | 0x40 | 0x20}\tctl\t1\t60\t50M\t=\t201\t250\t{REF[0:50]}\t{QUAL}\n"
    f"unflagged\t{1 | 0x80 | 0x10}\tctl\t201\t60\t50M\t=\t1\t-250\t{REF[200:250]}\t{QUAL}\n"
)
HEADER = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:ctl\tLN:300\n"


def rule_body(rule):
    text = SMK.read_text()
    start = text.index(f"rule {rule}:")
    end = text.find("\nrule ", start + 1)
    return text[start:end if end != -1 else None]


def methyldackel_command(rule):
    """Subcommand and long options of the rule's MethylDackel call."""
    m = re.search(r"MethylDackel (extract|mbias)(.*?)(?:\n\s*\"\"\"|$)", rule_body(rule), re.S)
    assert m, f"{rule} has no MethylDackel command"
    return m.group(1), re.findall(r"(?<!\S)--[A-Za-z]+", m.group(2))


@pytest.mark.parametrize("rule", RULES)
def test_rule_passes_keep_discordant(rule):
    _, options = methyldackel_command(rule)
    assert "--keepDiscordant" in options, f"{rule} drops pairs that lack the proper-pair flag"


@pytest.mark.skipif(
    shutil.which("MethylDackel") is None or shutil.which("samtools") is None,
    reason="MethylDackel or samtools not on PATH",
)
@pytest.mark.parametrize("rule", RULES)
def test_unflagged_pair_yields_calls(rule, tmp_path):
    fasta = tmp_path / "ctl.fa"
    fasta.write_text(f">ctl\n{REF}\n")
    subprocess.run(["samtools", "faidx", str(fasta)], check=True)
    sam = tmp_path / "pair.sam"
    sam.write_text(HEADER + PAIR)
    bam = tmp_path / "pair.bam"
    subprocess.run(["samtools", "view", "-b", "-o", str(bam), str(sam)], check=True)
    subprocess.run(["samtools", "index", str(bam)], check=True)

    subcommand, options = methyldackel_command(rule)
    if subcommand == "extract":
        prefix = tmp_path / "out"
        subprocess.run(
            ["MethylDackel", "extract", *options, str(fasta), str(bam), "-o", str(prefix)],
            check=True, capture_output=True,
        )
        rows = (tmp_path / "out_CpG.methylKit").read_text().splitlines()[1:]
    else:
        out = subprocess.run(
            ["MethylDackel", "mbias", *options, str(fasta), str(bam)],
            check=True, capture_output=True, text=True,
        ).stdout
        rows = out.splitlines()[1:]
    assert rows, f"{rule}: a pair without the proper-pair flag gave no methylation calls"
