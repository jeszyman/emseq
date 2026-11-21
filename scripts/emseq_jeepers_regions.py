#!/usr/bin/env python3
"""Write the CpG region file that JEEPERS's Bam2Txt.py reads, from a reference.

One row per CpG: contig, position minus 1, position plus 1, the 1-based
position of the C, a running index from 1, and a block number. JEEPERS
corrects each block in its own process. A new block starts where consecutive
CpGs are more than `gap` bases apart and the current block already holds
`min_block` CpGs.

--shift adds a constant to every position and --contig-name sets column 1,
to match a BAM written by emseq_jeepers_input.py with the same options.
"""
from __future__ import annotations

import argparse
import re

CG = re.compile("CG")


def cpg_rows(seq: str, contig: str, offset: int = 0, shift: int = 0, gap: int = 2000, min_block: int = 5000):
    """Yield (contig, pos - 1, pos + 1, pos, index, block) for each CG in seq.

    pos = offset + shift + i + 1 for a CG at string index i, where offset is
    the 0-based reference coordinate of the first base of seq.
    """
    block, in_block, last = 1, 0, None
    for index, m in enumerate(CG.finditer(seq.upper()), start=1):
        pos = offset + shift + m.start() + 1
        if last is not None and pos - last > gap and in_block >= min_block:
            block += 1
            in_block = 0
        in_block += 1
        last = pos
        yield (contig, pos - 1, pos + 1, pos, index, block)


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
    a = p.parse_args(argv)
    with pysam.FastaFile(a.fasta) as fa:
        contig, start, end = parse_region(a.region) if a.region else (fa.references[0], None, None)
        seq = fa.fetch(contig, start, end)
    n = block = 0
    with open(a.out, "w") as out:
        for row in cpg_rows(seq, a.contig_name or contig, offset=start or 0, shift=a.shift):
            out.write("\t".join(map(str, row)) + "\n")
            n, block = row[4], row[5]
    print(f"CpGs {n} blocks {block}")


if __name__ == "__main__":
    main()
