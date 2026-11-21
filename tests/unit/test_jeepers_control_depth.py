"""Depth per CpG from merge-table rows, for emseq_jeepers_control_depth.py."""
from emseq_jeepers_control_depth import cpg_depths, median_depth


def merge_line(positions):
    return "\t".join(["x"] * 10 + [positions] + ["x"] * 13) + "\n"


def test_depth_counts_rows_listing_each_cpg():
    depth = cpg_depths([merge_line("10,20"), merge_line("20,30"), merge_line("20"), merge_line("-")])
    assert depth == {"10": 1, "20": 3, "30": 1}


def test_median_is_over_covered_cpgs(tmp_path):
    table = tmp_path / "t_Merge.txt"
    table.write_text(merge_line("10,20,30") + merge_line("20,30") + merge_line("30"))
    assert median_depth(str(table)) == 2
