#!/usr/bin/env python3
"""Write the BAM that JEEPERS reads from a stored coordinate-sorted BAM.

JEEPERS (Wang 2026) accepts only read pairs flagged 99/147 or 163/83 whose
two reads align without gaps or clipping. This script keeps, within one
region, the pairs that are not flagged duplicate, have MAPQ 10 or higher on
both reads, carry one of those two flag pairs, and have a CIGAR of one M
operation on both reads. It writes a table of the pairs dropped by each
condition; a pair is counted under the first condition it fails, in the
order duplicate, mapq, flags, gapped. A read whose mate is not in the region,
and any secondary or supplementary record, is counted under flags.

For the spike-in controls three more options apply. --dedup-coordinates
keeps the first pair of each start, end, and strand. --shift adds a constant
to every position and to the contig length, so that all positions have the
same number of digits (JEEPERS sorts positions as text). --contig-name
renames the contig (JEEPERS needs the text "chr" in the name).

The output header holds only the region's contig.
"""
from __future__ import annotations

import argparse
import collections
import sys

MIN_MAPQ = 10
DUPLICATE = 0x400
# read 1 forward with read 2 reverse is an original-top pair; read 1 reverse
# with read 2 forward is an original-bottom pair
FLAG_PAIRS = {frozenset((99, 147)): "OT", frozenset((163, 83)): "OB"}
REASONS = ("duplicate", "mapq", "flags", "gapped", "coordinate_duplicate")
CIGAR_MATCH = 0


def ungapped(read) -> bool:
    cigar = read.cigartuples
    return cigar is not None and len(cigar) == 1 and cigar[0][0] == CIGAR_MATCH


def classify_pair(r1, r2) -> str:
    """Return "keep" or the first condition the pair fails."""
    if r1.is_duplicate or r2.is_duplicate:
        return "duplicate"
    if r1.mapping_quality < MIN_MAPQ or r2.mapping_quality < MIN_MAPQ:
        return "mapq"
    if frozenset((r1.flag, r2.flag)) not in FLAG_PAIRS:
        return "flags"
    if not (ungapped(r1) and ungapped(r2)):
        return "gapped"
    return "keep"


def pair_strand(r1, r2) -> str:
    return FLAG_PAIRS[frozenset((r1.flag, r2.flag))]


def select_pairs(reads, dedup_coordinates: bool = False):
    """Return (names to keep, Counter of pairs by outcome) for an iterable of reads.

    Reads are grouped by name in the order given. A name with other than two
    records is one pair dropped under "flags".
    """
    groups: dict[str, list] = collections.OrderedDict()
    for read in reads:
        groups.setdefault(read.query_name, []).append(read)
    keep: set[str] = set()
    counts: collections.Counter = collections.Counter()
    seen: set[tuple] = set()
    for name, members in groups.items():
        if len(members) != 2:
            if any(m.is_duplicate for m in members):
                counts["duplicate"] += 1
            else:
                counts["flags"] += 1
            continue
        reason = classify_pair(*members)
        if reason == "keep" and dedup_coordinates:
            key = (
                min(m.reference_start for m in members),
                max(m.reference_end for m in members),
                pair_strand(*members),
            )
            if key in seen:
                reason = "coordinate_duplicate"
            seen.add(key)
        counts[reason] += 1
        if reason == "keep":
            keep.add(name)
    return keep, counts


def output_header(header, contig: str, shift: int = 0, contig_name: str | None = None):
    """Header dict with only `contig`, its length raised by shift, renamed if asked."""
    d = header.to_dict()
    sq = [dict(s) for s in d.get("SQ", []) if s["SN"] == contig]
    if not sq:
        raise ValueError(f"contig {contig} is not in the BAM header")
    sq[0]["LN"] = int(sq[0]["LN"]) + shift
    if contig_name:
        sq[0]["SN"] = contig_name
    d["SQ"] = sq
    return d


def rewrite(read, header, shift: int = 0, contig_name: str | None = None):
    """The read against the output header, with POS and PNEXT raised by shift."""
    import pysam

    d = read.to_dict()
    if contig_name:
        d["ref_name"] = contig_name
    if d["next_ref_name"] not in ("=", "*"):
        d["next_ref_name"] = "="
    d["ref_pos"] = str(int(d["ref_pos"]) + shift)
    d["next_ref_pos"] = str(int(d["next_ref_pos"]) + shift)
    return pysam.AlignedSegment.from_dict(d, header)


def prepare(bam, region, out_bam, out_dropped, dedup_coordinates=False, shift=0, contig_name=None):
    import pysam

    with pysam.AlignmentFile(bam) as src:
        contig = region.split(":")[0] if region else src.references[0]
        fetch = (lambda: src.fetch(region=region)) if region else (lambda: src.fetch(contig))
        keep, counts = select_pairs(fetch(), dedup_coordinates)
        header = pysam.AlignmentHeader.from_dict(output_header(src.header, contig, shift, contig_name))
        with pysam.AlignmentFile(out_bam, "wb", header=header) as out:
            for read in fetch():
                if read.query_name in keep:
                    out.write(rewrite(read, header, shift, contig_name))
    pysam.index(out_bam)
    with open(out_dropped, "w") as fh:
        fh.write("outcome\tpairs\n")
        fh.write(f"keep\t{counts['keep']}\n")
        for reason in REASONS:
            fh.write(f"{reason}\t{counts[reason]}\n")
    return counts


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--bam", required=True)
    p.add_argument("--region", help="contig or contig:start-end; default the first contig")
    p.add_argument("--out-bam", required=True)
    p.add_argument("--out-dropped", required=True)
    p.add_argument("--dedup-coordinates", action="store_true")
    p.add_argument("--shift", type=int, default=0)
    p.add_argument("--contig-name")
    a = p.parse_args(argv)
    counts = prepare(a.bam, a.region, a.out_bam, a.out_dropped, a.dedup_coordinates, a.shift, a.contig_name)
    print("\t".join(f"{k}={v}" for k, v in sorted(counts.items())), file=sys.stderr)


if __name__ == "__main__":
    main()
