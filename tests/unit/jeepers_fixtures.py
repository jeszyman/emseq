"""Hand-built SAM records and BAM files shared by the JEEPERS test modules."""
import pysam

CONTIG = "ctg"
LENGTH = 5000
HEADER = f"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:{CONTIG}\tLN:{LENGTH}\n@SQ\tSN:other\tLN:100\n"


def sam_pair(name, start, end, strand="OT", read_len=50, mapq=(60, 60), dup=False, cigar=None, flags=None):
    """Two SAM lines for a fragment covering start..end (1-based, inclusive)."""
    left_flag, right_flag = flags or ((99, 147) if strand == "OT" else (163, 83))
    if dup:
        left_flag, right_flag = left_flag | 0x400, right_flag | 0x400
    read_len = min(read_len, end - start + 1)
    right = end - read_len + 1
    tlen = end - start + 1
    cig = cigar or f"{read_len}M"
    seq, qual = "A" * read_len, "I" * read_len
    return [
        (start, f"{name}\t{left_flag}\t{CONTIG}\t{start}\t{mapq[0]}\t{cig}\t=\t{right}\t{tlen}\t{seq}\t{qual}\n"),
        (right, f"{name}\t{right_flag}\t{CONTIG}\t{right}\t{mapq[1]}\t{read_len}M\t=\t{start}\t{-tlen}\t{seq}\t{qual}\n"),
    ]


def write_bam(path, pairs):
    """Write a coordinate-sorted, indexed BAM from a list of sam_pair() results."""
    lines = sorted((line for pair in pairs for line in pair), key=lambda x: x[0])
    sam = str(path) + ".sam"
    with open(sam, "w") as fh:
        fh.write(HEADER)
        fh.writelines(text for _, text in lines)
    with pysam.AlignmentFile(sam) as src, pysam.AlignmentFile(str(path), "wb", header=src.header) as out:
        for read in src:
            out.write(read)
    pysam.index(str(path))
    return str(path)


def reads_of(path):
    with pysam.AlignmentFile(str(path)) as bam:
        return list(bam.fetch(until_eof=True))
