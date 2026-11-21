# tests/unit/test_validation_end_handling.py
"""Known-answer tests for scripts/emseq_validation_end_handling.py.

The script turns MethylDackel mbias tables for the lambda and pUC19 controls
into the numbers the read-end validation report quotes. A position counts as
biased when the whole 95% Wilson interval of its pooled methylation lies more
than 5 percentage points from the pooled methylation of the middle 60% of
positions.
"""
import json
import subprocess
import sys

import pytest

from conftest import REPO_ROOT

import emseq_validation_end_handling as v

SCRIPT = REPO_ROOT / "scripts" / "emseq_validation_end_handling.py"
END_STRAND = {"r1_start": "OT", "r1_end": "OB", "r2_start": "OB", "r2_end": "OT"}


def test_wilson_interval():
    lo, hi = v.wilson(0, 10)
    assert lo == 0.0
    assert hi == pytest.approx(0.2775, abs=1e-4)
    lo, hi = v.wilson(50, 100)
    assert (lo, hi) == pytest.approx((0.4038, 0.5962), abs=1e-4)


def test_wilson_with_no_calls_is_none():
    assert v.wilson(0, 0) is None


def test_parse_bounds_reads_both_strands():
    line = "Suggested inclusion options: --OT 0,0,6,0 --OB 0,0,0,146"
    assert v.parse_bounds(line, "x.bounds.txt") == {"OT": (0, 0, 6, 0), "OB": (0, 0, 0, 146)}


def test_parse_bounds_keeps_a_strand_that_is_absent_as_none():
    assert v.parse_bounds("Suggested inclusion options: --OT 0,0,0,0", "x") == {"OT": (0, 0, 0, 0), "OB": None}


def test_parse_bounds_without_the_line_names_the_file():
    with pytest.raises(ValueError, match="lib1.puc19.CpG.bounds.txt"):
        v.parse_bounds("", "lib1.puc19.CpG.bounds.txt")


def test_end_bases_takes_each_end_from_the_strand_whose_position_1_is_that_end():
    bounds = {"OT": (3, 0, 7, 0), "OB": (2, 0, 6, 0)}
    assert v.end_bases(bounds, END_STRAND) == {"r1_start": 3, "r1_end": 2, "r2_start": 6, "r2_end": 7}


def test_end_bases_is_none_where_the_strand_has_no_bound():
    bounds = {"OT": (3, 0, 7, 0), "OB": None}
    assert v.end_bases(bounds, END_STRAND) == {"r1_start": 3, "r1_end": None, "r2_start": None, "r2_end": 7}


def test_reference_middle_pools_the_middle_sixty_percent():
    # 10 positions; the middle 60% are positions 3 to 8
    counts = [(0, 100)] * 2 + [(90, 10)] * 6 + [(0, 100)] * 2
    assert v.reference_middle(counts) == pytest.approx(0.90)


def test_biased_run_counts_consecutive_biased_positions_from_the_start():
    assert v.biased_run([(10, 90), (20, 80), (95, 5), (96, 4)], 0.95) == 2


def test_one_call_far_from_a_mid_range_reference_is_not_biased():
    # the interval of 0 of 1 is 0 to 0.79, which reaches within 5 points of 0.5
    assert v.biased_run([(0, 1), (50, 50)], 0.5) == 0


def test_one_unmethylated_call_against_a_high_reference_is_biased():
    # the same interval lies wholly below 0.95 - 0.05
    assert v.biased_run([(0, 1), (95, 5)], 0.95) == 1


def test_a_position_with_no_calls_ends_the_run():
    assert v.biased_run([(10, 90), (0, 0), (10, 90)], 0.95) == 1


def test_biased_within_counts_biased_positions_in_the_first_n():
    counts = [(95, 5), (10, 90), (95, 5), (10, 90)] + [(95, 5)] * 8
    assert v.biased_within(counts, 0.95, 10) == 2


def test_tlen_summary_gives_the_median_and_the_share_at_or_under_the_read_length(tmp_path):
    f = tmp_path / "x.tlen.tsv"
    # fragment lengths 100, 100, 150, 200: median 125, three of four at or under 150
    f.write_text("100\t2\n150\t1\n200\t1\n")
    assert v.tlen_summary([f], read_length=150) == {"median": 125.0, "share_le_read_pct": 75.0, "pairs": 4}


def test_tlen_summary_pools_libraries(tmp_path):
    a = tmp_path / "a.tlen.tsv"
    b = tmp_path / "b.tlen.tsv"
    a.write_text("100\t1\n")
    b.write_text("300\t3\n")
    assert v.tlen_summary([a, b], read_length=150) == {"median": 300.0, "share_le_read_pct": 25.0, "pairs": 4}


def series(start, n=150):
    """Counts per position: `start` at the first positions, 97% methylated after."""
    return list(start) + [(970, 30)] * (n - len(start))


def tables(puc19_r2_start, lambda_r1_start=()):
    low = [(1, 999)] * 150
    t = {}
    for context in ("CpG", "CHG", "CHH"):
        t[("unmeth_lambda", context)] = {end: list(low) for end in END_STRAND}
    t[("unmeth_lambda", "CpG")]["r1_start"] = list(lambda_r1_start) + low[len(lambda_r1_start):]
    t[("puc19", "CpG")] = {end: series(()) for end in END_STRAND}
    t[("puc19", "CpG")]["r2_start"] = series(puc19_r2_start)
    return t


def test_judge_fill_in_at_the_start_of_read_2_only():
    result = v.judge(tables(puc19_r2_start=[(100, 900)] * 4))
    assert result["P1"]["holds"] is True
    assert result["P1"]["biased_run"] == 4
    assert result["P2"]["holds"] is True
    assert result["P3"]["holds"] is True


def test_judge_p1_fails_without_a_biased_base_at_the_start_of_read_2():
    assert v.judge(tables(puc19_r2_start=()))["P1"]["holds"] is False


def test_judge_p2_and_p3_fail_when_lambda_is_biased_at_the_start_of_read_1():
    result = v.judge(tables(puc19_r2_start=[(100, 900)] * 4, lambda_r1_start=[(300, 700)] * 3))
    assert result["P2"]["holds"] is False
    assert result["P3"]["holds"] is False
    assert result["P3"]["biased"]["unmeth_lambda.CpG.r1_start"] == 3


MBIAS_HEADER = "Strand\tRead\tPosition\tnMethylated\tnUnmethylated\n"


def write_library(d, lib, source, context, rows, bounds):
    (d / f"{lib}.{source}.{context}.txt").write_text(MBIAS_HEADER + "".join(rows))
    (d / f"{lib}.{source}.{context}.bounds.txt").write_text(bounds)


def mbias_rows(r2_start_meth):
    rows = []
    for strand in ("OT", "OB"):
        for read in (1, 2):
            for pos in range(1, 21):
                meth = 970
                if (strand, read) == ("OB", 2) and pos <= 3:
                    meth = r2_start_meth
                rows.append(f"{strand}\t{read}\t{pos}\t{meth}\t{1000 - meth}\n")
    # position 21 has no row at all for OT read 1: zero calls there
    return rows


def run_script(tmp_path, bounds="Suggested inclusion options: --OT 0,0,0,0 --OB 0,0,4,0\n"):
    d = tmp_path / "mbias"
    d.mkdir()
    for lib in ("libA", "libB"):
        write_library(d, lib, "puc19", "CpG", mbias_rows(100), bounds)
        for context in ("CpG", "CHG", "CHH"):
            lam = [f"{s}\t{r}\t{p}\t1\t999\n" for s in ("OT", "OB") for r in (1, 2) for p in range(1, 21)]
            write_library(d, lib, "unmeth_lambda", context, lam, "Suggested inclusion options: --OT 0,0,0,0 --OB 0,0,0,0\n")
    c = tmp_path / "control"
    c.mkdir()
    for lib in ("libA", "libB"):
        for source, reads in (("puc19", 1200), ("unmeth_lambda", 50000)):
            (c / f"{lib}.{source}.reads.txt").write_text(f"{reads}\n")
            (c / f"{lib}.{source}.tlen.tsv").write_text("15\t3\n30\t1\n")
    out = tmp_path / "out"
    cmd = [sys.executable, str(SCRIPT), "--mbias-dir", str(d), "--control-dir", str(c), "--out-dir", str(out),
           "--library-index", "libA=2", "libB=1", "--methyldackel-version", "0.6.1", "--pipeline-commit", "abc1234"]
    return subprocess.run(cmd, capture_output=True, text=True), out


def test_script_writes_values_and_a_position_table_without_library_names(tmp_path):
    proc, out = run_script(tmp_path)
    assert proc.returncode == 0, proc.stderr
    values = json.loads((out / "values.json").read_text())
    assert values["n_libraries"] == 2
    assert values["methyldackel_version"] == "0.6.1"
    assert values["P1_holds"] is True
    assert values["puc19_r2_start_biased_run"] == 3
    assert values["puc19_r2_start_n_libraries_bound"] == 2
    assert values["puc19_r2_start_suggested_bound_max"] == 4
    assert values["puc19_r2_start_position1_pct"] == 10.0
    assert values["puc19_r2_start_position2_pct"] == 10.0
    # positions 4 to 20 of that end are 97% methylated
    assert values["puc19_r2_start_mean_after_run_pct"] == 97.0
    assert values["puc19_r1_start_mean_after_run_pct"] == 97.0
    assert values["puc19_usable_reads"] == 2400
    assert values["puc19_usable_reads_text"] == "2,400"
    assert values["unmeth_lambda_usable_reads_text"] == "100,000"
    assert values["read_length"] == 20
    assert values["puc19_tlen_median"] == 15.0
    assert values["puc19_share_tlen_le_read_pct"] == 75.0
    assert values["unmeth_lambda_overall_pct"] == 0.1
    assert values["unmeth_lambda_max_position_pct"] == 0.1
    assert values["puc19_total_calls_text"] == "160,000"
    table = (out / "position.tsv").read_text()
    assert "libA" not in table and "libB" not in table
    header = table.splitlines()[0].split("\t")
    assert header == ["source", "context", "strand", "read", "position", "library_index",
                      "n_methylated", "n_unmethylated", "lower", "upper"]
    assert {line.split("\t")[5] for line in table.splitlines()[1:]} == {"0", "1", "2"}


def test_script_stops_and_names_the_library_when_the_bounds_line_is_absent(tmp_path):
    proc, _ = run_script(tmp_path, bounds="")
    assert proc.returncode != 0
    assert "puc19.CpG.bounds.txt" in proc.stderr


# ---- plasma measures: loss by distance from the strand ends, calling settings

def test_end_loss_is_five_prime_minus_three_prime_with_its_interval():
    (loss, lower, upper), = v.end_loss([(900, 100)], [(700, 300)])
    assert loss == pytest.approx(0.20)
    assert (lower, upper) == pytest.approx((0.16605, 0.23395), abs=1e-4)


def test_end_loss_without_calls_at_one_end_is_none():
    assert v.end_loss([(0, 0)], [(700, 300)]) == [None]


def test_loss_run_counts_distances_whose_whole_interval_exceeds_the_threshold():
    losses = [(0.2, 0.17, 0.23), (0.1, 0.06, 0.14), (0.06, 0.04, 0.08), (0.2, 0.1, 0.3)]
    assert v.loss_run(losses) == 2


def test_mate_difference_pools_positions():
    difference, lower, upper = v.mate_difference([(450, 50), (450, 50)], [(880, 120)])
    assert difference == pytest.approx(0.02)
    assert (lower, upper) == pytest.approx((-0.00741, 0.04741), abs=1e-4)


def test_global_cpg_compares_settings_on_shared_cpgs(tmp_path):
    header = 'track type="bedGraph" description="x"\n'
    none = tmp_path / "none.bedGraph"
    neb = tmp_path / "neb.bedGraph"
    # CpG at 100 is called under both settings, CpG at 200 under "none" only
    none.write_text(header + "chr1\t100\t101\t50\t5\t5\nchr1\t200\t201\t0\t0\t10\n")
    neb.write_text(header + "chr1\t100\t101\t75\t6\t2\n")
    result = v.global_cpg({"none": none, "neb": neb})
    assert result["none"] == {"shared_pct": 50.0, "all_pct": 25.0, "calls": 20}
    assert result["neb"] == {"shared_pct": 75.0, "all_pct": 75.0, "calls": 8}


def test_paired_interval_of_per_library_differences():
    mean, lower, upper = v.paired_interval([1.0, 2.0, 3.0])
    assert mean == pytest.approx(2.0)
    # t(0.975, 2) = 4.303, standard error 1/sqrt(3)
    assert (lower, upper) == pytest.approx((-0.484, 4.484), abs=1e-3)


ENDS_HEADER = "fragments\tend\tmate\tdistance\tn_methylated\tn_unmethylated\tn_kept_neb\tn_kept_all10\n"


def write_ends(path, low_distances):
    rows = []
    for fragments in ("all", "long"):
        for end in ("5p", "3p"):
            for mate in (0, 1, 2):
                for distance in range(1, 81):
                    meth = 200 if end == "3p" and distance <= low_distances else 800
                    kept_neb = 500 if mate == 0 else 0
                    kept_all10 = 100 if mate == 0 else 0
                    rows.append(f"{fragments}\t{end}\t{mate}\t{distance}\t{meth}\t{1000 - meth}\t{kept_neb}\t{kept_all10}\n")
    path.write_text(ENDS_HEADER + "".join(rows))


def test_script_reports_the_plasma_loss_from_the_judged_libraries(tmp_path):
    proc, out = run_script(tmp_path)
    assert proc.returncode == 0, proc.stderr
    fragment_dir = tmp_path / "fragment"
    fragment_dir.mkdir()
    for lib in ("libA", "libB"):
        write_ends(fragment_dir / f"{lib}.host.ends.tsv", low_distances=3)
    cmd = [sys.executable, str(SCRIPT), "--mbias-dir", str(tmp_path / "mbias"), "--control-dir", str(tmp_path / "control"),
           "--out-dir", str(out), "--library-index", "libA=2", "libB=1", "--methyldackel-version", "0.6.1",
           "--pipeline-commit", "abc1234", "--fragment-dir", str(fragment_dir), "--fragment-sources", "host",
           "--check-index", "2"]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    values = json.loads((out / "values.json").read_text())
    assert values["host_n_libraries_judged"] == 1
    assert values["host_loss_run"] == 3
    assert values["host_loss_points_d1"] == 60.0
    assert values["host_loss_points_d10"] == 0.0
    assert values["host_check_passes"] is True
    assert values["host_check_difference_points"] == 0.0
    # read 1 and read 2 report the same methylation at the 3' end
    assert values["host_mate_difference_points"] == 0.0
    assert values["host_P5_holds"] is True
    # 500 of 1,000 once-per-fragment calls are kept under the NEB mask, 100 under the 10-base mask
    assert values["host_kept_neb_pct"] == 50.0
    assert values["host_kept_all10_pct"] == 10.0
    table = (out / "fragment_end.tsv").read_text().splitlines()
    assert table[0].split("\t")[:6] == ["source", "fragments", "end", "mate", "distance", "library_index"]
    assert {line.split("\t")[5] for line in table[1:]} == {"0", "1", "2"}
