library(argparse)
library(methylKit)

parser <- ArgumentParser()
parser$add_argument("--amp_file", required = TRUE)
parser$add_argument("--library_id", required = TRUE)
parser$add_argument("--treatment", type = "integer", required = TRUE)
parser$add_argument("--mincov", required = TRUE)
parser$add_argument("--build", required = TRUE)
parser$add_argument("--out_dir", required = TRUE)

args <- parser$parse_args()

obj <- methRead(
  location = args$amp_file,
  sample.id = args$library_id,
  assembly = args$build,
  treatment = args$treatment,
  pipeline = "amp",
  context = "CpG",
  resolution = "base",
  header = TRUE,
  mincov = args$mincov,
  dbtype = "tabix",
  dbdir = args$out_dir
)
