# tests/unit/test_sort_temp_prefix.py
"""Piped samtools sorts in one rule use distinct temporary-file prefixes.

samtools sort writes PREFIX.NNNN.bam when its input does not fit in memory.
Two sorts in one pipe that share a -T prefix fail with "File exists" once both
write temporary files, which happens at production depth but not on the test
libraries. emseq_dedup name-sorts before dupsifter and coordinate-sorts after
it, in one pipe.
"""
import random
import re
import shutil
import subprocess

import pytest

from conftest import REPO_ROOT

SMK = (REPO_ROOT / "workflows" / "emseq.smk").read_text()


def rule_bodies() -> dict:
    return dict(re.findall(r"^rule (\w+):\n(.*?)(?=^rule |\Z)", SMK, re.S | re.M))


def sort_prefixes(body: str) -> list:
    return re.findall(r'samtools sort[^\n]*?-T "([^"]+)"', body)


def test_each_rule_uses_distinct_sort_prefixes():
    for name, body in rule_bodies().items():
        prefixes = sort_prefixes(body)
        assert len(prefixes) == len(set(prefixes)), f"{name}: {prefixes}"


@pytest.mark.skipif(shutil.which("samtools") is None, reason="samtools not on PATH")
def test_dedup_sorts_with_rule_prefixes_spill_without_collision(tmp_path):
    prefixes = sort_prefixes(rule_bodies()["emseq_dedup"])
    assert len(prefixes) == 2
    name_tmp, coord_tmp = (p.replace("{params.temp_prefix}", str(tmp_path / "lib")) for p in prefixes)
    rng = random.Random(1)
    lines = ["@HD\tVN:1.6\tSO:unsorted", "@SQ\tSN:chr1\tLN:10000000"]
    seq, qual = "A" * 100, "I" * 100
    for i in range(200000):
        pos = rng.randint(1, 9000000)
        lines.append(f"r{i}\t99\tchr1\t{pos}\t60\t100M\t=\t{pos + 50}\t150\t{seq}\t{qual}")
        lines.append(f"r{i}\t147\tchr1\t{pos + 50}\t60\t100M\t=\t{pos}\t-150\t{seq}\t{qual}")
    sam = tmp_path / "in.sam"
    sam.write_text("\n".join(lines) + "\n")
    # -m 4M makes both sorts write temporary files, as at production depth
    cmd = (
        f'samtools view -bh "{sam}" '
        f'| samtools sort -n -m 4M -@ 2 -O BAM -T "{name_tmp}" -o - '
        f'| samtools sort -m 4M -@ 2 -O BAM -T "{coord_tmp}" -o "{tmp_path / "out.bam"}"'
    )
    run = subprocess.run(["bash", "-o", "pipefail", "-c", cmd], capture_output=True, text=True)
    assert "File exists" not in run.stderr
    assert run.returncode == 0, run.stderr[-2000:]
    count = subprocess.run(["samtools", "view", "-c", str(tmp_path / "out.bam")], capture_output=True, text=True)
    assert int(count.stdout) == 400000
