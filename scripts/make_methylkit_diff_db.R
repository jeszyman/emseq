#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
  library(methylKit)
})

# --- helpers ---
as_int <- function(x, nm) {
  xi <- suppressWarnings(as.integer(x))
  if (is.na(xi)) stop(sprintf("'%s' must be an integer (got: %s)", nm, x))
  xi
}
as_num_chunk <- function(x, nm="chunk_size") {
  xn <- suppressWarnings(as.numeric(x))
  if (is.na(xn)) stop(sprintf("'%s' must be numeric (accepts forms like 1e9 or 1000000000; got: %s)", nm, x))
  xn
}

# --- args ---
parser <- ArgumentParser()
parser$add_argument("--mbase", required=TRUE,
                    help="Path to methylBase_*.txt.bgz from unite step")
parser$add_argument("--cores", default="4",
                    help="Parallel workers (default: 4)")
parser$add_argument("--out_dir", required=TRUE,
                    help="Output directory for methylKit DBs")
parser$add_argument("--suffix", required=TRUE,
                    help="Suffix for DB files (e.g., 'test')")
parser$add_argument("--chunk_size", default="1e9",
                    help="Chunk size (accepts '1e9' style; default: 1e9)")

args <- parser$parse_args()

cores      <- as_int(args$cores, "cores")
chunk_size <- as_num_chunk(args$chunk_size, "chunk_size")

# --- load methylBase (tabix backend) ---

meth <- methylKit:::readMethylDB(args$mbase)

# --- Differential methylation ---
diff <- calculateDiffMeth(
  meth,
  mc.cores   = cores,
  chunk.size = chunk_size,
  save.db    = TRUE,
  dbdir      = args$out_dir
)

diff_path <- file.path(args$out_dir, sprintf("methylDiff_%s.txt.bgz", args$suffix))
message("Done. methylDiff: ", diff_path)
