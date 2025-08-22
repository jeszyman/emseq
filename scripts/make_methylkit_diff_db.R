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
  mincov = 10
)

# --- Unite ---
meth <- unite(merged_obj,
              destrand = FALSE,
              chunk.size = 1e9,
              mc.cores = as.numeric(args$cores),
              save.db = TRUE,
              min.per.group = 1,
              suffix = args$suffix,
              dbdir = args$out_dir)

# --- Diff methylation ---
diff <- calculateDiffMeth(meth,
                          mc.cores = as.numeric(args$cores),
                          chunk.size = 1e9,
                          save.db = TRUE,
                          dbdir = args$out_dir)

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

# --- Unite ---
meth <- unite(merged_obj,
              destrand = FALSE,
              chunk.size = 1e9,
              mc.cores = as.numeric(args$cores),
              save.db = TRUE,
              suffix = args$suffix,
              dbdir = args$out_dir)

# --- Diff methylation ---
diff <- calculateDiffMeth(meth,
                          mc.cores = as.numeric(args$cores),
                          chunk.size = 1e9,
                          save.db = TRUE,
                          dbdir = args$out_dir)

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

# --- Unite ---
meth <- unite(merged_obj,
              destrand = FALSE,
              chunk.size = 1e9,
              mc.cores = as.numeric(args$cores),
              save.db = TRUE,
              suffix = args$suffix,
              dbdir = args$out_dir)

# --- Diff methylation ---
diff <- calculateDiffMeth(meth,
                          mc.cores = as.numeric(args$cores),
                          chunk.size = 1e9,
                          save.db = TRUE,
                          dbdir = args$out_dir)
