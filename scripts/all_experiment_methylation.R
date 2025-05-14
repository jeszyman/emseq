library(argparse)
library(methylKit)
library(tidyverse)

# --- Argument Parsing ---
parser <- ArgumentParser()
parser$add_argument("--db_file", required = TRUE, help = "Path to tabix-indexed methylBase file")
parser$add_argument("--out_file", required = TRUE, help = "Output TSV path for percent methylation matrix")
parser$add_argument("--chunk_size", type = "double", default = 1e9, help = "Chunk size for methylKit operations")
args <- parser$parse_args()

# --- Check Header and Load Object ---
methylKit:::checkTabixHeader(args$db_file)
meth <- methylKit:::readMethylDB(args$db_file)

# --- Extract Percent Methylation Matrix ---
meth_matrix <- percMethylation(meth, rowids = TRUE, chunk.size = args$chunk_size)

# --- Write Output ---
write_tsv(as.data.frame(meth_matrix), args$out_file)
