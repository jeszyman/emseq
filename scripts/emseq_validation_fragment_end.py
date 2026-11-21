#!/usr/bin/env python3
"""Tally CpG calls by distance from the ends of the original strand.

End repair fills the 3' end of each original strand with unmethylated
cytosine, so methylation falls toward that end. The 5' end of the same strand
is not filled in and serves as the control. For each read pair this script
places each CpG call at its distance from the strand's 3' end and from its 5'
end, and writes the tallies by mate and once per fragment.

It also tallies, for CpGs covered by both mates of a pair, whether the mates
agree, and how many once-per-fragment calls two read-end masks leave in.

Exit codes: 0 success, 1 runtime error, 2 usage error.
"""
import argparse
import sys
from pathlib import Path

import pysam

# ==============================================================================
# Constants
# ==============================================================================

MIN_MAPQ = 10  # MethylDackel default
MIN_BASE_QUALITY = 5  # MethylDackel default
MAX_FRAGMENT = 1000  # as EMSEQ_PAIR_FILTER
# Bases at the start of read 2 that --nOT 0,0,0,5 --nOB 0,0,5,0 exclude, per
# strand (tests/unit/test_methyldackel_bounds.py)
NEB_MASK = {"OT": 4, "OB": 5}
CLIP = 10  # bases a 10-base clip removes at each read end
# On a fragment of this length or longer, a base within 80 of one strand end
# is nearer that end than the other, so the two ends can be compared
LONG_FRAGMENT = 160
SKIP_FLAGS = 0x4 | 0x8 | 0x100 | 0x200 | 0x400 | 0x800
FLUSH_EVERY = 200000  # reads between sweeps of unmatched mates


def parse_args():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument(
        "--bam", required=True, type=Path, help="coordinate-sorted, indexed alignment"
    )
    p.add_argument("--fasta", required=True, type=Path, help="indexed reference FASTA")
    p.add_argument("--out-ends", required=True, type=Path)
    p.add_argument("--out-agreement", required=True, type=Path)
    p.add_argument("--max-distance", type=int, default=150)
    p.add_argument(
        "--long-fragment",
        type=int,
        default=LONG_FRAGMENT,
        help="fragment length from which a fragment also enters the 'long' tally",
    )
    return p.parse_args()


def main():
    args = parse_args()
    try:
        every, agreement = tally(args.bam, args.fasta, max_distance=args.max_distance)
        long_only, _ = tally(
            args.bam,
            args.fasta,
            max_distance=args.max_distance,
            min_fragment=args.long_fragment,
        )
    except (OSError, ValueError) as err:
        print(f"error: {err}", file=sys.stderr)
        return 1
    write_ends({"all": every, "long": long_only}, args.out_ends)
    write_agreement(agreement, args.out_agreement)
    return 0


# ==============================================================================
# Tally
# ==============================================================================


def tally(bam, fasta, max_distance=150, min_fragment=0):
    """Tallies by strand end and by mate agreement for one alignment.

    Fragments shorter than min_fragment are left out.

    Returns (ends, agreement). ends maps (end, mate, distance) to a dict of
    n_methylated, n_unmethylated, n_kept_neb, and n_kept_all10, where end is
    "3p" or "5p", mate is 1, 2, or 0 for one call per fragment, and the two
    kept counts are filled for mate 0 only. agreement maps (mate, cycle) to a
    dict of n_agree and n_disagree.
    """
    ends, agreement = {}, {}
    pending = {}
    with pysam.AlignmentFile(str(bam)) as alignment, pysam.FastaFile(
        str(fasta)
    ) as reference:
        seen = 0
        for read in alignment:
            if not usable(read) or abs(read.template_length) < min_fragment:
                continue
            mate = pending.pop(read.query_name, None)
            calls = read_calls(read, reference)
            this = (read, calls)
            if mate is None:
                pending[read.query_name] = this
            else:
                add_fragment(ends, agreement, mate, this, max_distance)
            seen += 1
            if seen % FLUSH_EVERY == 0:
                sweep(pending, ends, agreement, read, max_distance)
        for single in pending.values():
            add_fragment(ends, agreement, single, None, max_distance)
    return ends, agreement


def usable(read):
    """Read of an inward-facing pair on one contig that MethylDackel would count."""
    if not read.is_paired or read.flag & SKIP_FLAGS:
        return False
    if read.mapping_quality < MIN_MAPQ or read.reference_id != read.next_reference_id:
        return False
    if read.is_reverse == read.mate_is_reverse:
        return False
    return 0 < abs(read.template_length) <= MAX_FRAGMENT


def strand_of(read):
    """Original top (OT) or original bottom (OB) strand the pair derives from."""
    return "OT" if read.is_read1 != read.is_reverse else "OB"


def read_calls(read, reference):
    """{reference position: (methylated, cycle, read length)} for the read's CpG calls.

    The position is that of the CpG cytosine on the top strand for an OT read
    and of the CpG guanine for an OB read. The cycle is the base's position in
    sequencing order, from 1.
    """
    strand = strand_of(read)
    start, end = read.reference_start, read.reference_end
    contig = read.reference_name
    left = max(0, start - 1)
    segment = reference.fetch(contig, left, end + 1).upper()
    offset = start - left
    if strand == "OT":
        sites = [
            start + i
            for i in range(end - start)
            if segment[offset + i] == "C"
            and segment[offset + i + 1 : offset + i + 2] == "G"
        ]
        methylated, unmethylated = "C", "T"
    else:
        sites = [
            start + i
            for i in range(end - start)
            if segment[offset + i] == "G"
            and offset + i > 0
            and segment[offset + i - 1] == "C"
        ]
        methylated, unmethylated = "G", "A"
    if not sites:
        return {}
    wanted = set(sites)
    sequence, qualities, length = (
        read.query_sequence,
        read.query_qualities,
        read.query_length,
    )
    calls = {}
    for query_pos, ref_pos in read.get_aligned_pairs(matches_only=True):
        if ref_pos not in wanted or qualities[query_pos] < MIN_BASE_QUALITY:
            continue
        base = sequence[query_pos]
        if base not in (methylated, unmethylated):
            continue
        cycle = length - query_pos if read.is_reverse else query_pos + 1
        calls[ref_pos] = (base == methylated, cycle, length)
    return calls


def fragment_bounds(read):
    """0-based, end-exclusive bounds of the fragment, from one mate."""
    length = abs(read.template_length)
    if read.is_reverse:
        return read.reference_end - length, read.reference_end
    return read.reference_start, read.reference_start + length


def add_fragment(ends, agreement, first, second, max_distance):
    """Add one fragment's calls; second is None when the mate was not seen."""
    reads = [first] + ([second] if second is not None else [])
    forward = next((r for r, _ in reads if not r.is_reverse), reads[0][0])
    start, end = fragment_bounds(forward)
    strand = strand_of(reads[0][0])
    by_mate = {1: {}, 2: {}}
    for read, calls in reads:
        by_mate[1 if read.is_read1 else 2] = calls
    for position in set(by_mate[1]) | set(by_mate[2]):
        from_left, from_right = position - start + 1, end - position
        distance = (
            {"5p": from_left, "3p": from_right}
            if strand == "OT"
            else {"5p": from_right, "3p": from_left}
        )
        one, two = by_mate[1].get(position), by_mate[2].get(position)
        pooled = one if one is not None else two
        # Read 1 supplies the call where both mates cover a base, so the NEB
        # mask removes a call only when read 2 alone covers it
        kept_neb = one is not None or two[1] > NEB_MASK[strand]
        kept_all10 = any(
            c is not None and CLIP < c[1] <= c[2] - CLIP for c in (one, two)
        )
        for end_name, d in distance.items():
            if not 1 <= d <= max_distance:
                continue
            for mate, call in ((1, one), (2, two)):
                if call is not None:
                    count(ends, (end_name, mate, d), call[0])
            row = count(ends, (end_name, 0, d), pooled[0])
            row["n_kept_neb"] += kept_neb
            row["n_kept_all10"] += kept_all10
        if one is not None and two is not None:
            agree = one[0] == two[0]
            for mate, call in ((1, one), (2, two)):
                cell = agreement.setdefault(
                    (mate, call[1]), {"n_agree": 0, "n_disagree": 0}
                )
                cell["n_agree" if agree else "n_disagree"] += 1


def count(ends, key, methylated):
    row = ends.setdefault(
        key,
        {"n_methylated": 0, "n_unmethylated": 0, "n_kept_neb": 0, "n_kept_all10": 0},
    )
    row["n_methylated" if methylated else "n_unmethylated"] += 1
    return row


def sweep(pending, ends, agreement, current, max_distance):
    """Tally reads whose mate can no longer arrive in a coordinate-sorted file."""
    stale = [
        name
        for name, (read, _) in pending.items()
        if read.reference_id != current.reference_id
        or read.reference_start < current.reference_start - MAX_FRAGMENT
    ]
    for name in stale:
        add_fragment(ends, agreement, pending.pop(name), None, max_distance)


# ==============================================================================
# Outputs
# ==============================================================================


def write_ends(tallies, path):
    """tallies maps a fragment set ("all" or "long") to its ends tally."""
    with open(path, "w") as out:
        out.write(
            "fragments\tend\tmate\tdistance\tn_methylated\tn_unmethylated\tn_kept_neb\tn_kept_all10\n"
        )
        for fragments, ends in tallies.items():
            for (end, mate, distance), row in sorted(ends.items()):
                out.write(
                    f"{fragments}\t{end}\t{mate}\t{distance}\t{row['n_methylated']}\t{row['n_unmethylated']}"
                    f"\t{row['n_kept_neb']}\t{row['n_kept_all10']}\n"
                )


def write_agreement(agreement, path):
    with open(path, "w") as out:
        out.write("mate\tcycle\tn_agree\tn_disagree\n")
        for (mate, cycle), cell in sorted(agreement.items()):
            out.write(f"{mate}\t{cycle}\t{cell['n_agree']}\t{cell['n_disagree']}\n")


if __name__ == "__main__":
    sys.exit(main())
