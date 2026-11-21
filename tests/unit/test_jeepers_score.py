"""Scoring of JEEPERS sibling-mode rows by emseq_jeepers_score.py, on hand-built rows."""
import pytest

import emseq_jeepers_score as sc


def row(start, end, strand, name, pos, before, after=None, marked="-", support="-", kind="Rd", contig="chr22"):
    """One 35-column sibling row; pos, before, and after are comma-separated text."""
    f = ["x"] * 35
    f[0], f[1], f[2], f[3], f[4], f[5] = contig, str(start), str(end), "1", strand, name
    f[9] = str(end - start + 1)
    f[10], f[11], f[24], f[25], f[26], f[34] = pos, before, after or before, marked, support, kind
    return "\t".join(f) + "\n"


def frag(*args, **kwargs):
    return sc.parse_row(row(*args, **kwargs))


def test_parse_row_reads_dash_marked_and_dash_support():
    f = frag(1001, 1200, "OT", "r", "1010,1190", "M,U")
    assert (f.contig, f.start, f.end, f.strand, f.name, f.kind) == ("chr22", 1001, 1200, "OT", "r", "Rd")
    assert f.pos == [1010, 1190] and f.before == ["M", "U"] and f.after == ["M", "U"]
    assert f.marked == [] and f.support == {}
    f = frag(1001, 1200, "OT", "r", "1010,1150,1190", "M,U,U", "M,M,U", "1150,1190", "3:4,-", "S")
    assert f.marked == [1150, 1190]
    assert f.support == {1150: 4, 1190: None}


def test_fallback_row_with_three_marked_cpgs_and_lone_dash_support():
    f = frag(1001, 1200, "OT", "r", "1150,1170,1190", "U,U,U", "m,m,q", "1150,1170,1190", "-", "Sf")
    assert f.marked == [1150, 1170, 1190]
    assert f.support == {1150: None, 1170: None, 1190: None}
    assert f.kind == "Sf"


def test_merge_row_of_24_columns_parses_with_after_equal_to_before():
    line = "\t".join(row(1001, 1200, "OB", "r", "1010", "U").rstrip("\n").split("\t")[:24])
    f = sc.parse_row(line)
    assert f.after == f.before == ["U"] and f.marked == [] and f.kind == ""


def test_ot_distance_end_minus_pos_plus_1():
    f = frag(1001, 1200, "OT", "r", "1200,1001", "M,M")
    assert sc.distances(f, 1200) == (1, 200)
    assert sc.distances(f, 1001) == (200, 1)
    assert sc.distances(f, 1191) == (10, 191)


def test_ob_distance_pos_plus_2_minus_start():
    # the called base is the G at pos + 1; with the C at 1000 the G is the fragment's first base
    f = frag(1001, 1200, "OB", "r", "1000,1199", "M,M")
    assert sc.distances(f, 1000) == (1, 200)
    assert sc.distances(f, 1199) == (200, 1)


def test_bins():
    assert [sc.bin_of(d) for d in (1, 10, 11, 30, 31, 50, 51, 80)] == [
        "1-10", "1-10", "11-20", "21-30", "31-50", "31-50", "51-80", "51-80"]
    assert sc.bin_of(81) is None and sc.bin_of(0) is None


def test_binary_drops_q_and_n():
    assert [sc.binary(s) for s in "MmUuqN"] == [1, 1, 0, 0, None, None]
    f = frag(1001, 1200, "OT", "r", "1100,1150,1190,1195", "M,U,U,U", "M,m,q,N", "1150,1190,1195", "-", "Sf")
    assert sc.method_calls(f, "jeepers") == {1100: 1, 1150: 1}
    acc = sc.new_acc()
    f = frag(1001, 1200, "OT", "r", "1010,1195,1198", "M,U,U", "M,q,N", "1195,1198", "-", "Sf")
    sc.add_gap(acc, f, methods=("jeepers",))
    assert sum(v[0] for k, v in acc.items() if k[0] == "jeepers") == 1


def test_majority_ignores_fragments_where_the_cpg_is_marked():
    fragments = [
        frag(1001, 1200, "OT", "target", "1190", "U", "M", "1190", "2:2", "S"),
        frag(1001, 1300, "OT", "a", "1190", "M"),
        frag(1001, 1300, "OT", "b", "1190", "M"),
        # three fragments read U at 1190 but have it marked, so they do not vote
        frag(1001, 1200, "OT", "c", "1190", "U", "U", "1190", "0:1", "S"),
        frag(1001, 1200, "OT", "d", "1190", "U", "U", "1190", "0:1", "S"),
        frag(1001, 1200, "OT", "e", "1190", "U", "U", "1190", "0:1", "S"),
    ]
    calls = sc.majority_calls(fragments)
    assert calls["target"] == {1190: "M"}
    assert set(calls) == {"target", "c", "d", "e"}
    tally = sc.tally_unmarked(fragments)
    assert sc.method_calls(fragments[0], "majority", tally) == {1190: 1}


def test_majority_tie_keeps_own_call():
    fragments = [
        frag(1001, 1200, "OT", "target", "1190", "U", "M", "1190", "1:2", "S"),
        frag(1001, 1300, "OT", "a", "1190", "M"),
        frag(1001, 1300, "OT", "b", "1190", "U"),
    ]
    assert sc.majority_calls(fragments)["target"] == {1190: "U"}
    # no other fragment at the CpG also keeps the own call
    assert sc.majority_calls(fragments[:1])["target"] == {1190: "U"}


def test_exclusion_drops_calls_within_30_of_3p_end():
    f = frag(1001, 1200, "OT", "r", "1010,1148,1149,1170,1171,1200", "M,M,M,U,U,U")
    assert set(sc.method_calls(f, "excl30")) == {1010, 1148, 1149, 1170}  # 1170 is 31 from the end
    assert set(sc.method_calls(f, "excl52")) == {1010, 1148}  # 1148 is 53 from the end
    b = frag(1001, 1200, "OB", "r", "1029,1030,1190", "U,U,M")
    assert set(sc.method_calls(b, "excl30")) == {1030, 1190}  # 1029 is 30 from the 3' end at the start


def test_drop_marked_removes_only_marked_calls_and_reports_partner_unmethylated_share():
    f = frag(1001, 1200, "OT", "r", "1010,1150,1195,1198", "M,M,U,U", "M,M,M,M", "1195,1198", "2:2,1:1", "S")
    d = sc.drop_marked(f)
    assert d.pos == [1010, 1150] and d.before == ["M", "M"] and d.marked == []
    assert sc.method_calls(f, "drop_marked") == {1010: 1, 1150: 1}
    acc = sc.new_acc()
    partner = {1010: 1, 1150: 1, 1195: 0, 1198: 1}
    sc.score_strand(acc, f, partner)
    # two calls dropped in the bin 1 to 10; the partner reads one of them unmethylated
    assert acc[("dropped_marked", "1-10")] == [2, 1]
    assert acc[("dropped_marked", "all")] == [2, 1]
    assert acc[("drop_marked", "all")] == [2, 2]
    assert acc[("uncorrected", "all")] == [4, 3]
    # both were changed to methylated: one on two reference fragments, one on one
    assert acc[("changed_to_methylated", "2")] == [1, 1]
    assert acc[("changed_to_methylated", "1")] == [1, 0]


def test_wilson_0_of_10_upper_about_27_8():
    est, lo, hi = sc.wilson(0, 10)
    assert est == 0 and lo == 0
    assert hi == pytest.approx(27.75, abs=0.05)


def test_newcombe_matches_a_worked_value():
    est, lo, hi = sc.newcombe(56, 70, 48, 80)
    assert est == pytest.approx(20.0, abs=0.01)
    assert lo == pytest.approx(5.2, abs=0.1)
    assert hi == pytest.approx(33.4, abs=0.1)


def test_wilson_of_0_calls_is_none_and_judges_unresolved():
    assert sc.wilson(0, 0) is None
    assert sc.newcombe(0, 0, 5, 10) is None and sc.newcombe(5, 10, 0, 0) is None
    for mode in ("at_most", "above", "at_least"):
        assert sc.judge(None, 0, mode) == "unresolved"


def test_judge_wide_interval_is_unresolved():
    assert sc.judge((50, 44, 53), 5, "at_most") == "unresolved"
    assert sc.judge((50, 47, 56), 5, "at_most") == "unresolved"
    assert sc.judge((50, 45, 55), 60, "at_most") == "met"  # exactly 5 points each side is resolved


def test_judge_at_most_met_when_lower_bound_at_or_below_limit():
    assert sc.judge((6, 5, 7), 5, "at_most") == "met"
    assert sc.judge((2, 1, 3), 5, "at_most") == "met"
    assert sc.judge((7, 5.1, 9), 5, "at_most") == "failed"


def test_judge_above_fails_when_interval_contains_0():
    assert sc.judge((1, -0.5, 2.5), 0, "above") == "failed"
    assert sc.judge((1, 0, 2), 0, "above") == "failed"
    assert sc.judge((2, 0.5, 3.5), 0, "above") == "met"


def test_judge_at_least_met_when_upper_bound_at_or_above_limit():
    assert sc.judge((-1, -2.5, 0), 0, "at_least") == "met"
    assert sc.judge((3, 1, 5), 0, "at_least") == "met"
    assert sc.judge((-3, -5, -0.1), 0, "at_least") == "failed"


def test_mid_value_uses_calls_over_30_from_both_ends():
    # 200-base OT fragment: 1030 is 30 from the 5' end, 1031 is 31; 1170 is 31 from the 3' end, 1171 is 30
    f = frag(1001, 1200, "OT", "r", "1030,1031,1100,1170,1171", "U,M,M,U,U", "M,M,M,M,M", "1170,1171", "-", "Sf")
    mid = sc.mid_value([f])
    assert mid["uncorrected"] == (2, 3)
    assert mid["jeepers"] == (3, 3)


def gap_fragments(n, meth_3p):
    """n OT fragments of 200 bases with one call 5 bases from each end; meth_3p of the 3' calls methylated."""
    return [
        frag(1001, 1200, "OT", f"r{i}", "1005,1196", "M," + ("M" if i < meth_3p else "U"))
        for i in range(n)
    ]


def test_gap_bin_met_when_interval_contains_5():
    acc = sc.new_acc()
    for f in gap_fragments(2000, 1890):  # 5' 100%, 3' 94.5%: gap 5.5 points, interval about 4.6 to 6.6
        sc.add_gap(acc, f, methods=("uncorrected",))
    r = sc.gap_rows(acc, ["uncorrected"])[0]
    assert (r["bin"], r["n_3p"], r["meth_3p"], r["n_5p"], r["meth_5p"]) == ("1-10", 2000, 1890, 2000, 2000)
    assert r["interval"][1] <= 5 <= r["interval"][2]
    assert r["verdict"] == "met"
    acc = sc.new_acc()
    for f in gap_fragments(2000, 1700):  # gap 15 points
        sc.add_gap(acc, f, methods=("uncorrected",))
    assert sc.gap_rows(acc, ["uncorrected"])[0]["verdict"] == "failed"
    # a bin with too few calls is unresolved, and so is a bin with none
    acc = sc.new_acc()
    for f in gap_fragments(20, 19):
        sc.add_gap(acc, f, methods=("uncorrected",))
    rows = sc.gap_rows(acc, ["uncorrected"])
    assert rows[0]["verdict"] == "unresolved"
    assert rows[1]["verdict"] == "unresolved" and rows[1]["interval"] is None


def test_gap_leaves_out_fragments_under_160_bases():
    acc = sc.new_acc()
    sc.add_gap(acc, frag(1001, 1159, "OT", "short", "1005,1155", "M,M"), methods=("uncorrected",))
    assert not acc
    sc.add_gap(acc, frag(1001, 1160, "OT", "long", "1005,1155", "M,M"), methods=("uncorrected",))
    assert acc[("uncorrected", "3p", "1-10")] == [1, 1] and acc[("uncorrected", "5p", "1-10")] == [1, 1]


def test_agreement_counts_only_calls_binary_in_both_strands():
    f = frag(1001, 1200, "OT", "r", "1150,1185,1195,1198", "M,U,U,U", "M,M,q,M", "1185,1195,1198", "2:2,-,1:1", "S")
    partner = {1150: 1, 1185: 1, 1195: 1}  # the partner has no call at 1198
    acc = sc.new_acc()
    sc.score_strand(acc, f, partner)
    assert acc[("uncorrected", "all")] == [3, 1]
    assert acc[("jeepers", "all")] == [2, 2]  # q at 1195 leaves the count
    assert acc[("uncorrected", "1-10")] == [1, 0] and acc[("jeepers", "1-10")] == [0, 0]
    assert acc[("uncorrected", "11-20")] == [1, 0] and acc[("jeepers", "11-20")] == [1, 1]
    assert acc[("uncorrected", "marked")] == [2, 0] and acc[("jeepers", "marked")] == [1, 1]
    assert acc[("excl30", "1-80")] == [1, 1] and acc[("uncorrected", "1-80")] == [3, 1]


def test_bar_fails_when_one_bin_fails():
    assert sc.bar_verdict(["met", "met", "failed", "unresolved"]) == "failed"
    assert sc.bar_verdict(["met", "met"]) == "met"
    assert sc.bar_verdict(["met", "unresolved"]) == "unresolved"
    assert sc.bar_verdict([]) == "unresolved"


def test_bar4_modes_by_bin():
    acc = sc.new_acc()
    for scope in ("1-10", "11-20", "21-30", "31-50"):
        acc[("jeepers", scope)] = [4000, 3600]
        acc[("uncorrected", scope)] = [4000, 3600]
    rows = {r["scope"]: r["verdict"] for r in sc.bar4_rows(acc)}
    # no difference: the near bins need an interval above 0, the far bins one that reaches 0
    assert rows == {"1-10": "failed", "11-20": "failed", "21-30": "met", "31-50": "met"}


def test_control_rows_compare_mid_fragment_with_bins_near_the_3p_end():
    fragments = [
        frag(10000001, 10000200, "OT", f"r{i}", "10000100,10000196", "M," + ("M" if i < 1900 else "U"))
        for i in range(2000)
    ]
    rows = {(r["method"], r["bin"]): r for r in sc.control_rows(fragments)}
    assert (rows[("uncorrected", "mid")]["n"], rows[("uncorrected", "mid")]["meth"]) == (2000, 2000)
    r = rows[("uncorrected", "1-10")]
    assert (r["n"], r["meth"]) == (2000, 1900)
    assert r["interval"][0] == pytest.approx(5.0)
    assert r["verdict"] == "met"
    assert rows[("uncorrected", "11-20")]["verdict"] == "unresolved"
