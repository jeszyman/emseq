"""The checks of emseq_jeepers_run.sh, driven with stand-in step scripts.

The stand-ins are four tiny scripts with the file names and options the
runner calls. They write fixed tables, which lets the runner's own checks
(empty table, lost rows, error text in a log) be exercised without JEEPERS.
"""
import os
import pathlib
import subprocess
import textwrap

RUNNER = pathlib.Path(__file__).resolve().parents[2] / "scripts" / "emseq_jeepers_run.sh"

ARGS = """
import argparse, sys
p = argparse.ArgumentParser()
for flag in ("--ibam", "--obam", "--ibed", "--ofile", "--cpu_num", "-c", "-d", "-i", "-o", "-n", "-m"):
    p.add_argument(flag)
a = p.parse_args()
"""

MERGE_ROWS = ["chr1\t1\t2\t1.0", "chr1\t3\t4\t1.0", "chr1\t5\t6\t2.0"]


def stub_dir(tmp_path, merge_rows=MERGE_ROWS, sibling_rows=MERGE_ROWS, correct_extra=""):
    d = tmp_path / "stub"
    (d / "FormatChange").mkdir(parents=True)
    (d / "FormatChange" / "bamFormatChange.py").write_text(ARGS + "open(a.obam, 'w').write('bam\\n')\n")
    (d / "Bam2Txt.py").write_text(ARGS + "open(a.ofile, 'w').write('r1\\nr2\\n')\n")
    (d / "MergeR1R2.py").write_text(
        ARGS + f"open(a.o, 'w').write(''.join(row + '\\n' for row in {merge_rows!r}))\n"
    )
    (d / "JEEPERS_correction.py").write_text(
        ARGS
        + textwrap.dedent(
            f"""
            open(a.o.replace('.txt', '_Sibling.txt'), 'w').write(''.join(row + '\\n' for row in {sibling_rows!r}))
            open(a.o.replace('.txt', '_args.txt'), 'w').write(str(a.m))
            {correct_extra}
            """
        )
    )
    return d


def run(tmp_path, setting="default", env_extra=None, **kwargs):
    d = stub_dir(tmp_path, **kwargs)
    (tmp_path / "in.bam").write_text("x")
    (tmp_path / "regions.bed").write_text("x")
    out = tmp_path / "out"
    env = dict(os.environ, SAMTOOLS="true", JEEPERS_PYTHON="python3", **(env_extra or {}))
    proc = subprocess.run(
        ["bash", str(RUNNER), str(tmp_path / "in.bam"), str(tmp_path / "regions.bed"), str(out), setting, str(d)],
        env=env,
        capture_output=True,
        text=True,
    )
    return proc, out


def test_runner_fails_on_empty_merge_table(tmp_path):
    proc, out = run(tmp_path, merge_rows=[])
    assert proc.returncode == 1
    assert "empty table" in proc.stderr
    assert not (out / "t_JEEPERS_Sibling.txt").exists()


def test_runner_fails_when_a_block_loses_rows(tmp_path):
    proc, out = run(tmp_path, sibling_rows=MERGE_ROWS[1:])
    assert proc.returncode == 1
    assert "differ" in proc.stderr
    assert (out / "blocks.tsv").read_text().split("\n")[1:3] == ["1.0\t2\t1", "2.0\t1\t1"]
    # a block missing from the output altogether also fails
    (tmp_path / "second").mkdir()
    proc, out = run(tmp_path / "second", sibling_rows=MERGE_ROWS[:2])
    assert proc.returncode == 1
    assert "2.0\t1\t0" in (out / "blocks.tsv").read_text()


def test_runner_fails_on_traceback_in_log(tmp_path):
    proc, _ = run(tmp_path, correct_extra="print('Traceback (most recent call last):')")
    assert proc.returncode == 1
    assert "Traceback or ERROR" in proc.stderr
    (tmp_path / "second").mkdir()
    proc, _ = run(tmp_path / "second", correct_extra="print('ERROR-1: Incomplete output for 1.0')")
    assert proc.returncode == 1


def test_runner_passes_and_writes_run_tsv(tmp_path):
    proc, out = run(tmp_path)
    assert proc.returncode == 0, proc.stderr
    rows = [line.split("\t") for line in (out / "run.tsv").read_text().rstrip("\n").split("\n")]
    assert rows[0] == ["step", "seconds", "rows"]
    assert [(r[0], r[2]) for r in rows[1:]] == [("format", "1"), ("bam2txt", "2"), ("merge", "3"), ("correct", "3")]
    assert (out / "t_Merge.txt").exists() and (out / "t_JEEPERS_Sibling.txt").exists()
    assert (out / "t_JEEPERS_args.txt").read_text() == "None"


def test_setting_0_9_passes_m_and_default_does_not(tmp_path):
    proc, out = run(tmp_path, setting="0.9")
    assert proc.returncode == 0, proc.stderr
    assert (out / "t_JEEPERS_args.txt").read_text() == "0.9"
    (tmp_path / "second").mkdir()
    proc, _ = run(tmp_path / "second", setting="0.5")
    assert proc.returncode == 2


def test_stop_after_merge_skips_the_correction(tmp_path):
    proc, out = run(tmp_path, env_extra={"JEEPERS_STOP_AFTER": "merge"})
    assert proc.returncode == 0, proc.stderr
    assert (out / "t_Merge.txt").exists() and not (out / "t_JEEPERS_Sibling.txt").exists()
