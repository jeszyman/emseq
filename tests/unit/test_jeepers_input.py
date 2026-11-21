"""Which read pairs emseq_jeepers_input.py passes to JEEPERS, on hand-built pairs."""
import pathlib
import sys

import pytest

pysam = pytest.importorskip("pysam")

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import emseq_jeepers_input as inp  # noqa: E402
from jeepers_fixtures import CONTIG, LENGTH, reads_of, sam_pair, write_bam  # noqa: E402


def classify(tmp_path, **kwargs):
    bam = write_bam(tmp_path / "one.bam", [sam_pair("p", 101, 300, **kwargs)])
    r1, r2 = reads_of(bam)
    return inp.classify_pair(r1, r2)


def test_pair_flagged_99_147_ungapped_is_kept(tmp_path):
    assert classify(tmp_path) == "keep"
    assert classify(tmp_path, strand="OB") == "keep"


def test_duplicate_flag_drops_pair(tmp_path):
    assert classify(tmp_path, dup=True) == "duplicate"


def test_mapq_9_on_one_read_drops_pair(tmp_path):
    assert classify(tmp_path, mapq=(60, 9)) == "mapq"


def test_flags_97_145_drop_pair(tmp_path):
    assert classify(tmp_path, flags=(97, 145)) == "flags"


def test_soft_clip_drops_pair(tmp_path):
    assert classify(tmp_path, cigar="5S45M") == "gapped"
    assert classify(tmp_path, cigar="20M1D30M") == "gapped"


def run(tmp_path, pairs, **kwargs):
    bam = write_bam(tmp_path / "in.bam", pairs)
    out, dropped = tmp_path / "out.bam", tmp_path / "dropped.tsv"
    inp.prepare(bam, CONTIG, str(out), str(dropped), **kwargs)
    table = dict(line.split("\t") for line in dropped.read_text().split("\n")[1:] if line)
    return reads_of(out), {k: int(v) for k, v in table.items()}, out


def test_dedup_coordinates_keeps_first_of_same_start_end_strand(tmp_path):
    pairs = [
        sam_pair("first", 101, 300),
        sam_pair("second", 101, 300),
        sam_pair("other_strand", 101, 300, strand="OB"),
        sam_pair("other_end", 101, 301),
    ]
    reads, table, _ = run(tmp_path, pairs, dedup_coordinates=True)
    assert {r.query_name for r in reads} == {"first", "other_strand", "other_end"}
    assert table["coordinate_duplicate"] == 1
    reads, table, _ = run(tmp_path, pairs)
    assert len(reads) == 8 and table["coordinate_duplicate"] == 0


def test_shift_raises_pos_pnext_and_contig_length(tmp_path):
    reads, _, out = run(tmp_path, [sam_pair("p", 5, 204)], shift=10_000_000)
    left, right = reads
    assert left.reference_start + 1 == 10_000_005
    assert left.next_reference_start + 1 == right.reference_start + 1 == 10_000_155
    assert right.next_reference_start + 1 == 10_000_005
    with pysam.AlignmentFile(str(out)) as bam:
        assert bam.lengths == (LENGTH + 10_000_000,)
    assert left.template_length == 200 and left.cigarstring == "50M"


def test_dropped_table_counts_each_reason(tmp_path):
    pairs = [
        sam_pair("keep1", 101, 300),
        sam_pair("keep2", 401, 600, strand="OB"),
        sam_pair("dup", 701, 900, dup=True),
        sam_pair("mapq", 1001, 1200, mapq=(5, 60)),
        sam_pair("flags", 1301, 1500, flags=(97, 145)),
        sam_pair("gapped", 1601, 1800, cigar="3S47M"),
        sam_pair("lone", 2001, 2200)[:1],
    ]
    reads, table, _ = run(tmp_path, pairs)
    assert table == {"keep": 2, "duplicate": 1, "mapq": 1, "flags": 2, "gapped": 1, "coordinate_duplicate": 0}
    assert sorted({r.query_name for r in reads}) == ["keep1", "keep2"]
    assert len(reads) == 4


def test_contig_name_renames_header_and_records(tmp_path):
    reads, _, out = run(tmp_path, [sam_pair("p", 101, 300)], contig_name="chrpUC19")
    with pysam.AlignmentFile(str(out)) as bam:
        assert bam.references == ("chrpUC19",)
    assert {r.reference_name for r in reads} == {"chrpUC19"}
    assert {r.next_reference_name for r in reads} == {"chrpUC19"}


def test_output_is_coordinate_sorted_and_indexed(tmp_path):
    reads, _, out = run(tmp_path, [sam_pair("b", 401, 600), sam_pair("a", 101, 300)])
    starts = [r.reference_start for r in reads]
    assert starts == sorted(starts)
    assert pathlib.Path(str(out) + ".bai").exists()
