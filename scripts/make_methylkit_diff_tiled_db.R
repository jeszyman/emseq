library(argparse)
library(methylKit)

# --- Argument Parsing ---
parser <- ArgumentParser()
parser$add_argument("--lib_db_list", required = TRUE)
parser$add_argument("--lib_id_list", required = TRUE)
parser$add_argument("--treatment_list", required = TRUE)
parser$add_argument("--cores", required = TRUE)
parser$add_argument("--out_dir", required = TRUE)
parser$add_argument("--suffix", required = TRUE)
parser$add_argument("--win_size", required = TRUE)
parser$add_argument("--step_size", required = TRUE)

args <- parser$parse_args()

lib_db_list <- unlist(strsplit(args$lib_db_list, " "))
lib_id_list <- unlist(strsplit(args$lib_id_list, " "))
treatment_list <- as.numeric(unlist(strsplit(args$treatment_list, " ")))

stopifnot(length(lib_db_list) == length(lib_id_list),
          length(lib_id_list) == length(treatment_list))

# --- Read methylation databases ---
merged_obj <- methRead(
  location = as.list(lib_db_list),
  sample.id = as.list(lib_id_list),
  treatment = treatment_list,
  context = "CpG",
  assembly = "hg38",
  dbtype = "tabix",
  mincov = 2
)


# --- Tile methylation ---
tiled_raw <- tileMethylCounts(
  merged_obj,
  win.size = as.numeric(args$win_size),
  step.size = as.numeric(args$step_size),
  cov.bases = 1,
  save.db = TRUE,
  suffix = args$suffix,
  dbdir = args$out_dir,
  sample.ids = lib_id_list,
  treatment = treatment_list,
  mc.cores = as.numeric(args$cores)
)

# --- Unite tiled windows into methylBaseDB ---
tiled_obj <- unite(
  tiled_raw,
  destrand = FALSE,
  save.db = TRUE,
  suffix = paste0(args$suffix, "_tiled"),
  dbdir = args$out_dir,
  mc.cores = as.numeric(args$cores)
)

# --- Diff methylation on tiles ---
diff <- calculateDiffMeth(
  tiled_obj,
  mc.cores = as.numeric(args$cores),
  chunk.size = 1e9,
  save.db = TRUE,
  dbdir = args$out_dir
)
