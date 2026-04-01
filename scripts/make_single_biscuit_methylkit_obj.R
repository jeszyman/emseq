library(argparse)
library(methylKit)

parser <- ArgumentParser()
parser$add_argument("--bismark_cov_bed", required = TRUE)
parser$add_argument("--library_id", required = TRUE)
parser$add_argument("--treatment", type = "integer", required = TRUE)
parser$add_argument("--out_dir", required = TRUE)

args <- parser$parse_args()

myobj= methRead(args$bismark_cov_bed,
                sample.id = args$library_id,
                treatment = 1,
                context="CpG",
                pipeline="bismarkCoverage",
                mincov = 2,
                assembly= "hg38",
                dbtype = "tabix",
                dbdir = args$out_dir)
