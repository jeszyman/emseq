#!/usr/bin/env python3
"""Hold out opposite strands of the same molecule before JEEPERS runs.

In a BAM written by emseq_jeepers_input.py, an original-top fragment and an
original-bottom fragment with identical start and end are the two strands of
one molecule, apart from the chance rate. This script lists those pairs,
puts one strand of each in set A and the other in set B by a seeded choice,
and writes the BAM without set A and the BAM without set B. A strand is then
scored from the JEEPERS run in which its partner was removed, so the
partner's calls cannot enter JEEPERS's reference for it.

It also lists near pairs, a bottom and a top fragment that share one end and
differ by one base at the other, as the chance comparison. Near pairs are
not removed from any BAM.

The pairs table has columns kind (same or near), top, bottom, and removed_in
(A or B: the set, and so the removed-from BAM, of the top read; empty for
near pairs).
"""
from __future__ import annotations

import argparse
import random

SEED = 20261120
TOP_FLAGS, BOTTOM_FLAGS = (99, 147), (163, 83)


def fragments(bam_path: str):
    """{(start, end, strand): read name} of the first fragment at each key, in file order.

    start and end are 1-based inclusive; strand is OT or OB.
    """
    import pysam

    spans: dict[str, list] = {}
    with pysam.AlignmentFile(bam_path) as bam:
        for read in bam.fetch(until_eof=True):
            strand = "OT" if read.flag in TOP_FLAGS else "OB" if read.flag in BOTTOM_FLAGS else None
            if strand is None:
                continue
            s = spans.setdefault(read.query_name, [read.reference_start + 1, read.reference_end, strand])
            s[0] = min(s[0], read.reference_start + 1)
            s[1] = max(s[1], read.reference_end)
    first: dict[tuple, str] = {}
    for name, (start, end, strand) in spans.items():
        first.setdefault((start, end, strand), name)
    return first


def strand_pairs(bam_path: str) -> list[tuple[str, str]]:
    """(top read name, bottom read name) for fragments with identical start and end."""
    first = fragments(bam_path)
    return [
        (name, first[(start, end, "OB")])
        for (start, end, strand), name in first.items()
        if strand == "OT" and (start, end, "OB") in first
    ]


def near_pairs(bam_path: str) -> list[tuple[str, str]]:
    """(top, bottom) where one end matches and the other differs by one base."""
    first = fragments(bam_path)
    out = []
    for (start, end, strand), name in first.items():
        if strand != "OT":
            continue
        for key in ((start, end + 1), (start, end - 1), (start + 1, end), (start - 1, end)):
            bottom = first.get((*key, "OB"))
            if bottom is not None:
                out.append((name, bottom))
    return out


def split(pairs, seed: int = SEED) -> tuple[set[str], set[str]]:
    """Put one read name of each pair in set A and the other in set B."""
    rng = random.Random(seed)
    a: set[str] = set()
    b: set[str] = set()
    for top, bottom in pairs:
        if rng.random() < 0.5:
            a.add(top), b.add(bottom)
        else:
            a.add(bottom), b.add(top)
    return a, b


def write_without(bam_path: str, out_path: str, names: set[str]) -> None:
    import pysam

    with pysam.AlignmentFile(bam_path) as src, pysam.AlignmentFile(out_path, "wb", header=src.header) as out:
        for read in src.fetch(until_eof=True):
            if read.query_name not in names:
                out.write(read)
    pysam.index(out_path)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--bam", required=True)
    p.add_argument("--out-a", required=True, help="BAM without the reads of set A")
    p.add_argument("--out-b", required=True, help="BAM without the reads of set B")
    p.add_argument("--out-pairs", required=True)
    p.add_argument("--seed", type=int, default=SEED)
    a = p.parse_args(argv)
    same, near = strand_pairs(a.bam), near_pairs(a.bam)
    set_a, set_b = split(same, a.seed)
    write_without(a.bam, a.out_a, set_a)
    write_without(a.bam, a.out_b, set_b)
    with open(a.out_pairs, "w") as fh:
        fh.write("kind\ttop\tbottom\tremoved_in\n")
        for top, bottom in same:
            fh.write(f"same\t{top}\t{bottom}\t{'A' if top in set_a else 'B'}\n")
        for top, bottom in near:
            fh.write(f"near\t{top}\t{bottom}\t\n")
    print(f"same {len(same)} near {len(near)}")


if __name__ == "__main__":
    main()
