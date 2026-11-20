# tests/unit/test_methyldackel_bounds.py
"""What MethylDackel 0.6.1 end bounds remove, on read pairs built by hand.

The reference has two contigs: ctgC is 300 C bases and ctgG is 300 G bases.
A pair on ctgC with read 1 forward is an original-top (OT) pair; a pair on
ctgG with read 1 reverse is an original-bottom (OB) pair. Read sequence equals
the reference, so every base gives one CHH call and a masked base is a
position with no call.

Facts pinned here:
- Bounds and mbias positions count along the stored alignment (reference
  left to right), not along the sequencing cycle. A reverse-strand read has
  its first sequenced base at its highest position.
- The first and third numbers of --nOT/--nOB mask that many bases at the
  alignment's left end. The second and fourth mask one base fewer at its
  right end.
- NEB's flags, --nOT 0,0,0,5 --nOB 0,0,5,0, both mask the first sequenced
  bases of read 2: 4 bases on OT pairs and 5 on OB pairs.
- A base masked in one mate is still called from the other mate where the
  mates overlap, so NEB's flags remove nothing on a fragment shorter than
  the read.
"""
import shutil
import subprocess

import pytest

pytestmark = pytest.mark.skipif(
    shutil.which("MethylDackel") is None or shutil.which("samtools") is None,
    reason="MethylDackel or samtools not on PATH",
)

# Flag strings for the three calling settings compared in the read-end
# validation. all10 masks 10 bases at each end of both reads on both strands;
# the right-end numbers are 11 because they mask one base fewer than given.
END_SETTINGS = {
    "none": "",
    "neb": "--nOT 0,0,0,5 --nOB 0,0,5,0",
    "all10": "--nOT 10,11,10,11 --nOB 10,11,10,11",
}

# Strand whose mbias table position 1 sits at each read end. The table counts
# from the alignment's left end, which is the first sequenced base of a
# forward read and the last sequenced base of a reverse read.
END_STRAND = {"r1_start": "OT", "r1_end": "OB", "r2_start": "OB", "r2_end": "OT"}

HEADER = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:ctgC\tLN:300\n@SQ\tSN:ctgG\tLN:300\n"
BASE = {"ctgC": "C", "ctgG": "G"}
# 0x1 paired, 0x40/0x80 read 1/2, 0x10 reverse, 0x20 mate reverse
OT_R1, OT_R2 = 1 | 0x40 | 0x20, 1 | 0x80 | 0x10
OB_R1, OB_R2 = 1 | 0x40 | 0x10, 1 | 0x80 | 0x20


def record(flag, contig, pos, length, mate_pos, tlen, seq=None):
    seq = seq or BASE[contig] * length
    return f"p\t{flag}\t{contig}\t{pos}\t60\t{length}M\t=\t{mate_pos}\t{tlen}\t{seq}\t{'I' * length}\n"


PAIRS = {
    # 250 bp fragments, mates apart
    "ot_apart": record(OT_R1, "ctgC", 1, 50, 201, 250) + record(OT_R2, "ctgC", 201, 50, 1, -250),
    "ob_apart": record(OB_R2, "ctgG", 1, 50, 201, 250) + record(OB_R1, "ctgG", 201, 50, 1, -250),
    # 80 bp fragment, mates overlap at 31-50
    "ot_overlap": record(OT_R1, "ctgC", 1, 50, 31, 80) + record(OT_R2, "ctgC", 31, 50, 1, -80),
    # 50 bp fragments read through by both mates
    "ot_full": record(OT_R1, "ctgC", 1, 50, 1, 50) + record(OT_R2, "ctgC", 1, 50, 1, -50),
    "ob_full": record(OB_R2, "ctgG", 1, 50, 1, 50) + record(OB_R1, "ctgG", 1, 50, 1, -50),
    # reads of 30 bases, as after adapter trimming
    "ot_short": record(OT_R1, "ctgC", 1, 30, 201, 230) + record(OT_R2, "ctgC", 201, 30, 1, -230),
    # one unmethylated base at the left end of each alignment
    "ot_mark": record(OT_R1, "ctgC", 1, 50, 201, 250, "T" + "C" * 49)
    + record(OT_R2, "ctgC", 201, 50, 1, -250, "T" + "C" * 49),
    "ob_mark": record(OB_R2, "ctgG", 1, 50, 201, 250, "A" + "G" * 49)
    + record(OB_R1, "ctgG", 201, 50, 1, -250, "A" + "G" * 49),
}


@pytest.fixture(scope="module")
def workdir(tmp_path_factory):
    d = tmp_path_factory.mktemp("bounds")
    (d / "ref.fa").write_text(">ctgC\n" + "C" * 300 + "\n>ctgG\n" + "G" * 300 + "\n")
    subprocess.run(["samtools", "faidx", str(d / "ref.fa")], check=True)
    for name, records in PAIRS.items():
        (d / f"{name}.sam").write_text(HEADER + records)
        subprocess.run(["samtools", "view", "-b", "-o", str(d / f"{name}.bam"), str(d / f"{name}.sam")], check=True)
        subprocess.run(["samtools", "index", str(d / f"{name}.bam")], check=True)
    return d


def called_positions(workdir, flags, pair):
    """1-based positions with a CHH call, mapped to their call depth."""
    subprocess.run(
        ["MethylDackel", "extract", "--CHH", "--keepDiscordant", *flags.split(),
         str(workdir / "ref.fa"), str(workdir / f"{pair}.bam"), "-o", str(workdir / "out")],
        check=True, capture_output=True,
    )
    depth = {}
    for line in (workdir / "out_CHH.bedGraph").read_text().splitlines()[1:]:
        f = line.split("\t")
        depth[int(f[1]) + 1] = int(f[4]) + int(f[5])
    return depth


def removed(workdir, flags, pair):
    return sorted(set(called_positions(workdir, "", pair)) - set(called_positions(workdir, flags, pair)))


def mbias_rows(workdir, pair):
    """(strand, read, position, n_methylated, n_unmethylated) per table row."""
    out = subprocess.run(
        ["MethylDackel", "mbias", "--noSVG", "--CHH", "--noCpG", "--keepDiscordant",
         str(workdir / "ref.fa"), str(workdir / f"{pair}.bam")],
        check=True, capture_output=True, text=True,
    ).stdout
    rows = []
    for line in out.splitlines()[1:]:
        strand, read, pos, meth, unmeth = line.split("\t")
        rows.append((strand, int(read), int(pos), int(meth), int(unmeth)))
    return rows


def span(a, b):
    return list(range(a, b + 1))


def test_version_is_the_pinned_one():
    out = subprocess.run(["MethylDackel", "--version"], capture_output=True, text=True)
    assert (out.stdout + out.stderr).split()[0] == "0.6.1"


@pytest.mark.parametrize("flags, pair, expected", [
    # left-end numbers mask that many bases
    ("--nOT 5,0,0,0", "ot_apart", span(1, 5)),
    ("--nOT 0,0,5,0", "ot_apart", span(201, 205)),
    ("--nOB 0,0,5,0", "ob_apart", span(1, 5)),
    ("--nOB 5,0,0,0", "ob_apart", span(201, 205)),
    # right-end numbers mask one base fewer
    ("--nOT 0,5,0,0", "ot_apart", span(47, 50)),
    ("--nOT 0,0,0,5", "ot_apart", span(247, 250)),
    ("--nOB 0,5,0,0", "ob_apart", span(247, 250)),
    ("--nOB 0,0,0,5", "ob_apart", span(47, 50)),
    # the same ends are masked on reads of 30 bases
    ("--nOT 0,0,0,5", "ot_short", span(227, 230)),
    ("--nOT 5,0,0,0", "ot_short", span(1, 5)),
])
def test_which_positions_a_bound_removes(workdir, flags, pair, expected):
    assert removed(workdir, flags, pair) == expected


def test_neb_flags_mask_the_start_of_read_2_four_bases_on_ot_and_five_on_ob(workdir):
    # OT read 2 is the reverse read at 201-250; its first sequenced base is 250.
    assert removed(workdir, END_SETTINGS["neb"], "ot_apart") == span(247, 250)
    # OB read 2 is the forward read at 1-50; its first sequenced base is 1.
    assert removed(workdir, END_SETTINGS["neb"], "ob_apart") == span(1, 5)


@pytest.mark.parametrize("pair", ["ot_apart", "ob_apart"])
def test_all10_masks_ten_bases_at_each_of_the_four_ends(workdir, pair):
    assert removed(workdir, END_SETTINGS["all10"], pair) == span(1, 10) + span(41, 50) + span(201, 210) + span(241, 250)


def test_overlapped_bases_are_called_once(workdir):
    assert set(called_positions(workdir, "", "ot_overlap").values()) == {1}


def test_a_base_masked_in_read_2_is_called_from_read_1_where_mates_overlap(workdir):
    # read 2 spans 31-80; masking its left 5 bases (31-35) removes nothing,
    # because read 1 (1-50) covers them
    assert removed(workdir, "--nOT 0,0,5,0", "ot_overlap") == []


@pytest.mark.parametrize("pair", ["ot_full", "ob_full"])
def test_neb_flags_remove_nothing_on_a_fragment_shorter_than_the_read(workdir, pair):
    assert removed(workdir, END_SETTINGS["neb"], pair) == []


@pytest.mark.parametrize("pair, strand", [("ot_mark", "OT"), ("ob_mark", "OB")])
def test_mbias_position_counts_from_the_left_end_of_the_alignment(workdir, pair, strand):
    unmethylated = [(s, read, pos) for s, read, pos, meth, unmeth in mbias_rows(workdir, pair) if unmeth]
    # the marked base is the leftmost aligned base of both mates: position 1,
    # also for the reverse-strand mate, whose leftmost base was sequenced last
    assert unmethylated == [(strand, 1, 1), (strand, 2, 1)]


def test_end_strand_names_the_strand_whose_position_1_is_that_read_end():
    # forward reads: OT read 1 and OB read 2 (position 1 = first sequenced base)
    # reverse reads: OB read 1 and OT read 2 (position 1 = last sequenced base)
    assert END_STRAND == {"r1_start": "OT", "r1_end": "OB", "r2_start": "OB", "r2_end": "OT"}
