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
as_num <- function(x, nm) {
  xn <- suppressWarnings(as.numeric(x))
  if (is.na(xn)) stop(sprintf("'%s' must be numeric (e.g., 1e9 or 1000000000; got: %s)", nm, x))
  xn
}

# --- args ---
parser <- ArgumentParser()
parser$add_argument("--lib_db_list", required = TRUE)
parser$add_argument("--lib_id_list", required = TRUE)
parser$add_argument("--treatment_list", required = TRUE)
parser$add_argument("--cores", required = TRUE)
parser$add_argument("--out_dir", required = TRUE)
parser$add_argument("--suffix", required = TRUE)
parser$add_argument("--assembly", required = TRUE)
parser$add_argument("--mincov", required = TRUE)
parser$add_argument("--win_size", required = TRUE)
parser$add_argument("--chunk_size", required = TRUE)
parser$add_argument("--min_per_group", default="1", help="Min samples per group per CpG (integer; default: 1)")
args <- parser$parse_args()

dir.create(args$out_dir, recursive = TRUE, showWarnings = FALSE)

# --- parse ---
lib_db_list    <- split_ws(args$lib_db_list)
lib_id_list    <- split_ws(args$lib_id_list)
treatment_list <- as.numeric(split_ws(args$treatment_list))

if (length(lib_db_list) != length(lib_id_list) ||
    length(lib_id_list) != length(treatment_list)) {
  stop(sprintf("Length mismatch: lib_db_list=%d, lib_id_list=%d, treatment_list=%d",
               length(lib_db_list), length(lib_id_list), length(treatment_list)))
}

cores      <- as_int(args$cores, "cores")
mincov     <- as_int(args$mincov, "mincov")
win_size   <- as_int(args$win_size, "win_size")      # FIX: was args$winsize → numeric(0)
min_per_grp  <- as_int(args$min_per_group, "min_per_group")
chunk_size <- as_num(args$chunk_size, "chunk_size")

# --- destructive cleanup: only *_tiled* (both base & diff) ---
tiled_patterns <- c(
  file.path(args$out_dir, sprintf("methylBase_%s_tiled*.txt.bgz*", args$suffix)),
  file.path(args$out_dir, sprintf("methylDiff_%s_tiled*.txt.bgz*", args$suffix))
)
tiled_paths <- unlist(lapply(tiled_patterns, Sys.glob))
if (length(tiled_paths) > 0) suppressWarnings(file.remove(tiled_paths))

# --- read methylation DBs ---
merged_obj <- methRead(
  location  = as.list(lib_db_list),
  sample.id = as.list(lib_id_list),
  treatment = treatment_list,
  context   = "CpG",
  assembly  = args$assembly,
  dbtype    = "tabix",
  mincov    = mincov
)

# --- tile per sample ---
tiled_raw <- tileMethylCounts(
  merged_obj,
  win.size   = win_size,
  step.size  = win_size,
  cov.bases  = 1,
  save.db    = TRUE,
  suffix     = args$suffix,
  dbdir      = args$out_dir,
  mc.cores   = cores
)

# --- unite tiles into methylBaseDB ---
tiled_obj <- unite(
  tiled_raw,
  destrand = FALSE,
  save.db  = TRUE,
  suffix   = paste0(args$suffix, "_tiled"),
  dbdir    = args$out_dir,
  min.per.group = min_per_grp,
  mc.cores = cores
)

# --- differential methylation on tiles ---
diff <- calculateDiffMeth(
  tiled_obj,
  mc.cores   = cores,
  chunk.size = chunk_size,   # FIX: missing comma previously
  save.db    = TRUE,
  dbdir      = args$out_dir,
)
