library(argparse)
library(methylKit)
library(tidyverse)
library(Rsamtools)
library(GenomicRanges)

# --- Argument Parsing ---
parser <- ArgumentParser()
parser$add_argument("--db_file", required = TRUE, help = "Path to tabix-indexed methylBase file")
parser$add_argument("--out_file", required = TRUE, help = "Output TSV path for percent methylation matrix")
parser$add_argument("--chunk_size", type = "double", default = 1e9, help = "Chunk size for methylKit operations")
args <- parser$parse_args()

# --- Check Header and Load Object ---
methylKit:::checkTabixHeader(args$db_file)
meth <- methylKit:::readMethylDB(args$db_file)

chroms <- Rsamtools::seqnamesTabix(args$db_file)

first_chunk <- TRUE

for (chrom in chroms) {
    message(sprintf("Processing chromosome: %s", chrom))

    chrom_gr <- GRanges(seqnames = chrom, ranges = IRanges(start = 1, end = .Machine$integer.max))

    meth_sub <- selectByOverlap(meth, chrom_gr)

    if (is.null(meth_sub) || nrow(meth_sub) == 0) next

    mat <- percMethylation(meth_sub, rowids = FALSE, chunk.size = args$chunk_size)
    df <- as.data.frame(mat)

    chunk_data <- getData(meth_sub)
    df$coord <- paste(chunk_data$chr, chunk_data$start, chunk_data$end, sep = ".")

    df <- df %>% select(coord, everything())

    write.table(df, file = args$out_file, append = !first_chunk,
                col.names = first_chunk, row.names = FALSE, sep = "\t", quote = FALSE)
    first_chunk <- FALSE
}
