#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
  library(methylKit)
})

# --- helpers ---
split_ws <- function(x) trimws(unlist(strsplit(x, "\\s+")))
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
parser$add_argument("--lib_db_list", required=TRUE, help="Space-separated tabix files")
parser$add_argument("--lib_id_list", required=TRUE, help="Space-separated sample IDs")
parser$add_argument("--treatment_list", required=TRUE, help="Space-separated 0/1 indicators")
parser$add_argument("--cores", default="4", help="Parallel workers (default: 4)")
parser$add_argument("--out_dir", required=TRUE, help="Output directory for methylKit DBs")
parser$add_argument("--suffix", required=TRUE, help="Suffix for DB files (e.g., 'test')")
parser$add_argument("--assembly", default="hg38", help="Genome assembly (default: hg38)")
parser$add_argument("--mincov", default="10", help="Minimum coverage per CpG (integer; default: 10)")
parser$add_argument("--min_per_group", default="1", help="Min samples per group per CpG (integer; default: 1)")
parser$add_argument("--chunk_size", default="1e9", help="Chunk size (accepts '1e9' style; default: 1e9)")

args <- parser$parse_args()

# --- parse & validate ---
lib_db_list    <- split_ws(args$lib_db_list)
lib_id_list    <- split_ws(args$lib_id_list)
treatment_list <- as.numeric(split_ws(args$treatment_list))

if (length(lib_db_list) != length(lib_id_list) ||
    length(lib_id_list) != length(treatment_list)) {
  stop(sprintf("Length mismatch: lib_db_list=%d, lib_id_list=%d, treatment_list=%d",
               length(lib_db_list), length(lib_id_list), length(treatment_list)))
}

cores        <- as_int(args$cores, "cores")
mincov       <- as_int(args$mincov, "mincov")
min_per_grp  <- as_int(args$min_per_group, "min_per_group")
chunk_size   <- as_num_chunk(args$chunk_size, "chunk_size")

# destructive overwrite of methylBase only
mbase_path <- file.path(args$out_dir, sprintf("methylBase_%s.txt.bgz", args$suffix))
suppressWarnings(file.remove(mbase_path, paste0(mbase_path, ".tbi")))

# --- Read (tabix + CpG only) ---
merged_obj <- methRead(
  location   = as.list(lib_db_list),
  sample.id  = as.list(lib_id_list),
  treatment  = treatment_list,
  context    = "CpG",
  assembly   = args$assembly,
  dbtype     = "tabix",
  mincov     = mincov
)

# --- Unite (destrand hardcoded TRUE, save.db always TRUE) ---
meth <- unite(
  merged_obj,
  destrand      = TRUE,
  chunk.size    = chunk_size,
  mc.cores      = cores,
  save.db       = TRUE,
  min.per.group = min_per_grp,
  suffix        = args$suffix,
  dbdir         = args$out_dir
)

message("Done. methylBase: ", mbase_path)
