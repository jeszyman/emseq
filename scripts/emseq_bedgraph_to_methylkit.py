#!/usr/bin/env python3
"""Convert a MethylDackel CpG bedGraph to the table MethylDackel writes with --methylKit.

The bedGraph holds one line per cytosine of a CpG (chrom, 0-based start, end,
percent, methylated, unmethylated) without the strand. The strand is read from
the reference: the call sits on a C for the top strand (F) and on a G for the
bottom strand (R). Also prints, on stderr, the autosomal totals of methylated
and unmethylated calls and the number of positions.

Usage: emseq_bedgraph_to_methylkit.py <in.bedGraph[.gz]> <ref.fa> <out.methylKit>
"""
from __future__ import annotations

import gzip
import sys

AUTOSOMES = {str(i) for i in range(1, 23)} | {f"chr{i}" for i in range(1, 23)}


def read_fai(path: str) -> dict:
    """{contig: (length, offset, bases per line, bytes per line)} from a .fai index."""
    index = {}
    with open(path) as fh:
        for line in fh:
            name, length, offset, linebases, linewidth = line.split("\t")[:5]
            index[name] = (int(length), int(offset), int(linebases), int(linewidth))
    return index


def contig_sequence(fasta, entry) -> bytes:
    """Upper-case sequence of one contig, read through its index entry."""
    length, offset, linebases, linewidth = entry
    lines = -(-length // linebases)
    fasta.seek(offset)
    raw = fasta.read(length + lines * (linewidth - linebases))
    return raw.replace(b"\n", b"").replace(b"\r", b"").upper()[:length]


def strand_of(base: int) -> str | None:
    """F for a call on C, R for a call on G, None for another base."""
    if base == 67:
        return "F"
    if base == 71:
        return "R"
    return None


def convert(bedgraph, fasta_path: str, out) -> dict:
    index = read_fai(fasta_path + ".fai")
    totals = {"positions": 0, "meth": 0, "unmeth": 0, "other_base": 0}
    current, seq = None, b""
    out.write("chrBase\tchr\tbase\tstrand\tcoverage\tfreqC\tfreqT\n")
    with open(fasta_path, "rb") as fasta:
        for line in bedgraph:
            if line.startswith("track"):
                continue
            chrom, start, _, _, meth, unmeth = line.rstrip("\n").split("\t")
            if chrom != current:
                current, seq = chrom, contig_sequence(fasta, index[chrom])
            start, meth, unmeth = int(start), int(meth), int(unmeth)
            strand = strand_of(seq[start])
            if strand is None:
                totals["other_base"] += 1
                continue
            cov = meth + unmeth
            pos = start + 1
            out.write(f"{chrom}.{pos}\t{chrom}\t{pos}\t{strand}\t{cov}\t{100 * meth / cov:.2f}\t{100 * unmeth / cov:.2f}\n")
            if chrom in AUTOSOMES:
                totals["positions"] += 1
                totals["meth"] += meth
                totals["unmeth"] += unmeth
    return totals


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 3:
        print(__doc__, file=sys.stderr)
        return 2
    opener = gzip.open if argv[0].endswith(".gz") else open
    with opener(argv[0], "rt") as bedgraph, open(argv[2], "w") as out:
        totals = convert(bedgraph, argv[1], out)
    print("\t".join(f"{k}={v}" for k, v in totals.items()), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
