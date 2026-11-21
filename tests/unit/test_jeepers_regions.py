"""The CpG region rows emseq_jeepers_regions.py writes, on short sequences."""
from emseq_jeepers_regions import cpg_rows, parse_region


def test_cg_at_index_0_gives_position_1():
    assert list(cpg_rows("CGA", "c")) == [("c", 0, 2, 1, 1, 1)]


def test_offset_and_shift_add_to_position():
    rows = list(cpg_rows("AACG", "c", offset=100, shift=10_000_000))
    assert rows == [("c", 10_000_102, 10_000_104, 10_000_103, 1, 1)]


def test_new_block_needs_gap_over_2000_and_5000_cpgs():
    # CpGs at 1, 3, 20, 31, 42: the gap of 17 after one block of 2 opens
    # block 2; the gap of 11 after a block of 1 does not; the gap of 11 after
    # a block of 2 does.
    seq = "CGCG" + "A" * 15 + "CG" + "A" * 9 + "CG" + "A" * 9 + "CG"
    rows = list(cpg_rows(seq, "c", gap=10, min_block=2))
    assert [r[3] for r in rows] == [1, 3, 20, 31, 42]
    assert [r[5] for r in rows] == [1, 1, 2, 2, 3]
    # a gap of exactly 10 is not over 10
    seq = "CGCG" + "A" * 7 + "CG"
    assert [r[5] for r in cpg_rows(seq, "c", gap=10, min_block=2)] == [1, 1, 1]
    # defaults: no split in a short sequence
    assert {r[5] for r in cpg_rows("CG" + "A" * 3000 + "CG", "c")} == {1}


def test_index_runs_from_1():
    assert [r[4] for r in cpg_rows("CGACGTTCG", "c")] == [1, 2, 3]


def test_lowercase_sequence_is_read():
    assert [r[3] for r in cpg_rows("acgtCgcG", "c")] == [2, 5, 7]


def test_parse_region():
    assert parse_region("chr22") == ("chr22", None, None)
    assert parse_region("chr21:20,000,001-25000000") == ("chr21", 20_000_000, 25_000_000)


def spaced(gaps):
    """A sequence with a CG, then one more CG after each gap of that many bases between C and C."""
    return "CG" + "".join("A" * (g - 2) + "CG" for g in gaps)


def test_max_block_caps_block_size_and_cuts_at_the_widest_gap_near_the_cap():
    # 10 CpGs; gaps between consecutive CpGs: 5 5 5 9 5 5 5 5 5
    seq = spaced([5, 5, 5, 9, 5, 5, 5, 5, 5])
    rows = list(cpg_rows(seq, "c", max_block=6, cut_window=4))
    # the first block would take a 7th CpG; of the gaps before CpGs 4 to 7 the widest (9) is before CpG 5
    assert [r[5] for r in rows] == [1, 1, 1, 1, 2, 2, 2, 2, 2, 2]
    assert [r[4] for r in rows] == list(range(1, 11))
    # with the wide gap outside the window, the cut falls on the latest of the equal gaps, at the cap
    rows = list(cpg_rows(seq, "c", max_block=6, cut_window=2))
    assert [r[5] for r in rows] == [1, 1, 1, 1, 1, 1, 2, 2, 2, 2]


def test_max_block_is_never_exceeded_and_default_is_unchanged():
    seq = spaced([3] * 99)
    blocks = [r[5] for r in cpg_rows(seq, "c", max_block=7, cut_window=3)]
    sizes = [blocks.count(b) for b in sorted(set(blocks))]
    assert max(sizes) <= 7 and sum(sizes) == 100
    assert blocks == sorted(blocks)
    assert {r[5] for r in cpg_rows(seq, "c")} == {1}
    # the gap rule still opens a block before the cap is reached
    seq = spaced([3, 3, 20, 3, 3, 3])
    assert [r[5] for r in cpg_rows(seq, "c", gap=10, min_block=2, max_block=5)] == [1, 1, 1, 2, 2, 2, 2]
