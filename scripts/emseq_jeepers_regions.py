#!/usr/bin/env python3
"""Write the CpG region file that JEEPERS's Bam2Txt.py reads, from a reference.

One row per CpG: contig, position minus 1, position plus 1, the 1-based
position of the C, a running index from 1, and a block number. JEEPERS
corrects each block in its own process. A new block starts where consecutive
CpGs are more than `gap` bases apart and the current block already holds
`min_block` CpGs. --max-block also caps the number of CpGs in a block, cutting
at the widest gap near the cap; JEEPERS holds one array of fragments by
CpGs per block, so its memory grows with the square of block size.

--shift adds a constant to every position and --contig-name sets column 1,
to match a BAM written by emseq_jeepers_input.py with the same options.
"""
from __future__ import annotations

import argparse
import re

CG = re.compile("CG")


def cpg_rows(seq: str, contig: str, offset: int = 0, shift: int = 0, gap: int = 2000, min_block: int = 5000,
             max_block: int | None = None, cut_window: int = 2000):
    """Yield (contig, pos - 1, pos + 1, pos, index, block) for each CG in seq.

    pos = offset + shift + i + 1 for a CG at string index i, where offset is
    the 0-based reference coordinate of the first base of seq.

    With max_block, no block holds more than max_block CpGs: when a block
    would take one more, it is cut at the widest gap between consecutive
    CpGs among its last cut_window CpGs (the later gap on a tie), and the
    CpGs after the cut open the next block. JEEPERS gives a fragment to
    each block in which it covers a CpG, so cutting at a wide gap keeps
    fragments from being split between two blocks.
    """
    positions = [offset + shift + m.start() + 1 for m in CG.finditer(seq.upper())]
    starts = [0] if positions else []  # index of the first CpG of each block
    for i in range(1, len(positions)):
        held = i - starts[-1]
        if positions[i] - positions[i - 1] > gap and held >= min_block:
            starts.append(i)
        elif max_block is not None and held >= max_block:
            first = max(starts[-1] + 1, i - cut_window + 1)
            cut = max(range(first, i + 1), key=lambda c: (positions[c] - positions[c - 1], c))
            starts.append(cut)
    block, next_start = 0, 0
    for i, pos in enumerate(positions):
        if next_start < len(starts) and starts[next_start] == i:
            block += 1
            next_start += 1
        yield (contig, pos - 1, pos + 1, pos, i + 1, block)


def parse_region(region: str):
    """(contig, 0-based start or None, end or None) from contig or contig:start-end."""
    if ":" not in region:
        return region, None, None
    contig, span = region.rsplit(":", 1)
    start, end = span.replace(",", "").split("-")
    return contig, int(start) - 1, int(end)


def main(argv=None):
    import pysam

    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--fasta", required=True)
    p.add_argument("--region", help="contig or contig:start-end; default the first contig")
    p.add_argument("--out", required=True)
    p.add_argument("--shift", type=int, default=0)
    p.add_argument("--contig-name")
    p.add_argument("--max-block", type=int, help="largest number of CpGs in a block; default no limit")
    a = p.parse_args(argv)
    with pysam.FastaFile(a.fasta) as fa:
        contig, start, end = parse_region(a.region) if a.region else (fa.references[0], None, None)
        seq = fa.fetch(contig, start, end)
    n = block = 0
    with open(a.out, "w") as out:
        for row in cpg_rows(seq, a.contig_name or contig, offset=start or 0, shift=a.shift, max_block=a.max_block):
            out.write("\t".join(map(str, row)) + "\n")
            n, block = row[4], row[5]
    print(f"CpGs {n} blocks {block}")


if __name__ == "__main__":
    main()
