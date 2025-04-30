library(argparse)
library(methylKit)

parser <- ArgumentParser()
parser$add_argument("--bismark_cov_bed", required = TRUE)
parser$add_argument("--library_id", required = TRUE)
parser$add_argument("--treatment", type = "integer", required = TRUE)
parser$add_argument("--mincov", type = "integer", required = TRUE)
parser$add_argument("--out_dir", required = TRUE)
parser$add_argument("--build", required = TRUE)

args <- parser$parse_args()

myobj= methRead(args$bismark_cov_bed,
                sample.id = args$library_id,
                treatment = args$treatment,
                context="CpG",
                pipeline="bismarkCoverage",
                mincov = args$mincov,
                assembly=args$build,
                dbtype = "tabix",
                dbdir = args$out_dir)
