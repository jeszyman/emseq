# tests/unit/test_spike_alignment.py
"""Control reads are measured in an alignment to the sample reference and the controls together.

Aligned to a control alone, sample reads match it partially and are counted as
control methylation. In 10.1 million read pairs of one public sheared-DNA library,
pUC19 read 82.8% methylated when aligned alone and 96.9% when aligned with the
sample reference.
The shell of each rule in workflows/emseq.smk is run here with a stub aligner.
"""
import ast
import gzip
import os
import re
import shutil
import subprocess
from types import SimpleNamespace as NS

import pytest

from conftest import REPO_ROOT

SMK = (REPO_ROOT / "workflows" / "emseq.smk").read_text()
SQ = "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:100000\n@SQ\tSN:lam\tLN:48502\n@SQ\tSN:puc\tLN:2686\n"
SEQ, QUAL = "A" * 50, "I" * 50

needs_samtools = pytest.mark.skipif(shutil.which("samtools") is None, reason="samtools not on PATH")


def rule_body(name):
    match = re.search(rf"^rule {name}:\n(.*?)(?=^rule |^def |\Z)", SMK, re.S | re.M)
    assert match, f"rule {name} not found"
    return match.group(1)


def rule_shell(name):
    match = re.search(r'shell:\s*\n\s*(""".*?""")', rule_body(name), re.S)
    assert match, f"rule {name} has no shell block"
    # Snakemake reads the block as a Python string, so its escapes are resolved
    return ast.literal_eval(match.group(1))


def run_shell(name, tmp_path, path_prefix=None, **fields):
    """Run a rule's shell in bash strict mode, as Snakemake does."""
    script = rule_shell(name).format(log=NS(cmd=tmp_path / "log"), threads=1, **fields)
    env = dict(os.environ)
    if path_prefix:
        env["PATH"] = f"{path_prefix}:{env['PATH']}"
    result = subprocess.run(["bash", "-c", f"set -euo pipefail\n{script}"], env=env, capture_output=True, text=True)
    log = (tmp_path / "log").read_text() if (tmp_path / "log").exists() else ""
    return result.returncode, result.stderr + log


def record(name, flag, rname, pos, rnext="=", pnext=None, mapq=60):
    return f"{name}\t{flag}\t{rname}\t{pos}\t{mapq}\t50M\t{rnext}\t{pnext or pos}\t0\t{SEQ}\t{QUAL}\n"


# 0x1 paired, 0x4 unmapped, 0x8 mate unmapped, 0x10/0x20 reverse, 0x40/0x80 read 1/2, 0x800 supplementary
ALIGNED = (
    record("human", 0x61, "chr1", 100) + record("human", 0x91, "chr1", 200)
    + record("lambda", 0x61, "lam", 100) + record("lambda", 0x91, "lam", 200)
    + record("puc", 0x61, "puc", 100) + record("puc", 0x91, "puc", 200)
    + record("split_pair", 0x61, "lam", 300, "chr1", 300) + record("split_pair", 0x91, "chr1", 300, "lam", 300)
    + record("half", 0x49, "lam", 400) + f"half\t{0x85}\tlam\t400\t0\t*\t=\t400\t0\t{SEQ}\t{QUAL}\n"
    + record("lambda", 0x861, "lam", 500)
)


@pytest.fixture
def stub_bin(tmp_path):
    """A bwameth.py that prints the SAM of ALIGNED for an alignment and touches nothing for an index."""
    d = tmp_path / "bin"
    d.mkdir()
    sam = tmp_path / "aligned.sam"
    sam.write_text(SQ + ALIGNED)
    stub = d / "bwameth.py"
    stub.write_text(f'#!/bin/bash\nif [[ "$1" == index* ]]; then exit 0; fi\ncat "{sam}"\n')
    stub.chmod(0o755)
    return d


def fai(path, *contigs):
    path.write_text("".join(f"{c}\t1000\t0\t60\t61\n" for c in contigs))
    return path


def names(bam):
    out = subprocess.run(["samtools", "view", bam], capture_output=True, text=True, check=True).stdout
    return [(l.split("\t")[0], l.split("\t")[2]) for l in out.splitlines()]


def align_controls(tmp_path, stub_bin):
    """Run the align rule on ALIGNED; return the BAM of every control."""
    fais = [fai(tmp_path / "lam.fai", "lam"), fai(tmp_path / "puc.fai", "puc")]
    out = tmp_path / "controls.bam"
    code, err = run_shell(
        "emseq_align_bwameth_spike", tmp_path, stub_bin,
        wildcards=NS(library_id="lib"), params=NS(temp_prefix=tmp_path / "tmp" / "lib"),
        input=NS(ref="x.fa", r1="r1", r2="r2", fais=" ".join(map(str, fais))), output=NS(bam=out),
    )
    assert code == 0, err
    return out


@needs_samtools
def test_alignment_keeps_primary_pairs_on_a_control(tmp_path, stub_bin):
    assert sorted(names(align_controls(tmp_path, stub_bin))) == [
        ("half", "lam"), ("lambda", "lam"), ("lambda", "lam"), ("puc", "puc"), ("puc", "puc")]


@needs_samtools
@pytest.mark.parametrize("control,contig,kept", [
    ("lambda", "lam", [("half", "lam"), ("lambda", "lam"), ("lambda", "lam")]),
    ("puc19", "puc", [("puc", "puc"), ("puc", "puc")]),
])
def test_split_writes_exactly_the_reads_of_one_control(tmp_path, stub_bin, control, contig, kept):
    out = tmp_path / "one.bam"
    code, err = run_shell(
        "emseq_split_spike", tmp_path, stub_bin,
        wildcards=NS(library_id="lib", emseq_ref_name=control),
        input=NS(bam=align_controls(tmp_path, stub_bin), fai=fai(tmp_path / "c.fai", contig)), output=NS(bam=out),
    )
    assert code == 0, err
    assert sorted(names(out)) == kept
    header = subprocess.run(["samtools", "view", "-H", out], capture_output=True, text=True, check=True).stdout
    assert re.findall(r"^@SQ\tSN:(\S+)", header, re.M) == [contig]


def run_index(tmp_path, stub_bin, fastas):
    target = tmp_path / "combined.fa"
    return target, run_shell(
        "emseq_bwa_meth_index", tmp_path, stub_bin,
        wildcards=NS(emseq_ref_name="x"), params=NS(fasta_target=target),
        input=NS(fastas=" ".join(map(str, fastas))),
    )


@needs_samtools
def test_index_joins_gzipped_and_plain_fastas_without_a_final_newline(tmp_path, stub_bin):
    host = tmp_path / "host.fa.gz"
    with gzip.open(host, "wt") as handle:
        handle.write(">chr1 x\nACGT\nACGT")
    control = tmp_path / "lam.fa"
    control.write_text(">lam y\nCCGG\n")
    target, (code, err) = run_index(tmp_path, stub_bin, [host, control])
    assert code == 0, err
    assert target.read_text() == ">chr1 x\nACGT\nACGT\n>lam y\nCCGG\n"
    assert [l.split("\t")[0] for l in (tmp_path / "combined.fa.fai").read_text().splitlines()] == ["chr1", "lam"]


@needs_samtools
def test_index_stops_when_a_control_repeats_a_contig_name(tmp_path, stub_bin):
    for name in ("a.fa", "b.fa"):
        (tmp_path / name).write_text(">chr1\nACGT\n")
    _, (code, err) = run_index(tmp_path, stub_bin, [tmp_path / "a.fa", tmp_path / "b.fa"])
    assert code != 0
    assert "contig name" in err


def index_fastas(name):
    """Input FASTAs of the index of one reference, from the function of the workflow."""
    scope = dict(spike_host_ref="host", spike_builds=["lam", "puc"], D_INPUTS="in",
                 EMSEQ_REF_INPUTS={"host": "host.fa.gz", "lam": "lam.fa", "puc": "puc.fa"})
    exec(re.search(r"^EMSEQ_SPIKE_REF = .*$", SMK, re.M).group(0), scope)
    exec(re.search(r"^def _emseq_index_fastas\(.*?(?=^rule )", SMK, re.S | re.M).group(0), scope)
    return scope["_emseq_index_fastas"](NS(emseq_ref_name=name))


def test_combined_index_holds_the_sample_reference_and_every_control():
    assert index_fastas("host_spike") == ["in/host.fa.gz", "in/lam.fa", "in/puc.fa"]


def test_index_of_one_reference_holds_that_reference():
    assert index_fastas("lam") == ["in/lam.fa"]


def test_the_align_rule_reads_the_combined_reference():
    body = rule_body("emseq_align_bwameth_spike")
    assert "{EMSEQ_SPIKE_REF}/{EMSEQ_SPIKE_REF}.fa" in body and "{{emseq_ref_name}}" not in body


@pytest.mark.parametrize("wrapper", ["test.smk", "test-analysis.smk", "validation.smk"])
def test_every_wrapper_names_the_sample_reference_of_the_controls(wrapper):
    text = (REPO_ROOT / "workflows" / wrapper).read_text()
    assert re.search(r"^spike_host_ref\s*=", text, re.M)
    assert text.index("spike_host_ref") < text.index('include: "emseq.smk"')
