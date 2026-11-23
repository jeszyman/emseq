#!/usr/bin/env Rscript
# Compare the three read-end settings on the full-depth cohort.
#
# Usage: emseq_validation_settings_compare.R <work dir> <out dir> <index map>
#
# Reads, from the work dir:
#   global/<library>.<setting>.txt   autosomal totals written by emseq_bedgraph_to_methylkit.py
#   dmr/<setting>/methylDiff_*.txt.bgz   tile-level differential methylation, patients against healthy
#   index map: TSV with columns library and library_index (leading '#' lines
#   ignored), so the tables name libraries only by their index
# Writes, to the out dir:
#   global.tsv      one row per library index and setting: positions, methylated, unmethylated, percent
#   global_diff.tsv per setting pair: paired difference in percent over libraries, with 95% interval,
#                   and per setting: patient minus healthy difference with 95% interval
#   tiles.tsv       per setting: tiles tested and tiles called at two thresholds
#   tiles_pair.tsv  per setting pair: shared tiles, correlation and largest change of the methylation
#                   difference, calls shared, Jaccard index, calls whose sign differs
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
work <- args[1]; out <- args[2]; index_map <- args[3]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
settings <- c("none", "neb", "all10")
pairs <- list(c("neb", "none"), c("all10", "none"), c("all10", "neb"))

# ---- global methylation ----
read_global <- function(path) {
  kv <- strsplit(strsplit(readLines(path, n = 1), "\t")[[1]], "=")
  vals <- setNames(as.numeric(sapply(kv, `[`, 2)), sapply(kv, `[`, 1))
  parts <- strsplit(sub("\\.txt$", "", basename(path)), ".", fixed = TRUE)[[1]]
  data.table(library = parts[1], setting = parts[2], positions = vals[["positions"]],
             meth = vals[["meth"]], unmeth = vals[["unmeth"]])
}
glob <- rbindlist(lapply(Sys.glob(file.path(work, "global", "*.txt")), read_global))
glob[, pct := 100 * meth / (meth + unmeth)]
glob[, group := ifelse(grepl("^NH_", library), "healthy", "patient")]
idx <- fread(index_map, skip = "library\tlibrary_index")
missing <- setdiff(unique(glob$library), idx$library)
if (length(missing)) stop("index map has no entry for: ", paste(missing, collapse = ", "))
glob[idx, library_index := i.library_index, on = "library"]
fwrite(glob[order(library_index, setting), .(library_index, setting, group, positions, meth, unmeth, pct)],
       file.path(out, "global.tsv"), sep = "\t")

wide <- dcast(glob, library + group ~ setting, value.var = "pct")
ci <- function(x) { t <- t.test(x); c(mean(x), t$conf.int[1], t$conf.int[2]) }
rows <- list()
for (p in pairs) {
  d <- wide[[p[1]]] - wide[[p[2]]]
  v <- ci(d)
  rows[[length(rows) + 1]] <- data.table(comparison = paste(p[1], "minus", p[2]), n = length(d),
                                         estimate = v[1], lower = v[2], upper = v[3])
}
for (s in settings) {
  t <- t.test(wide[group == "patient"][[s]], wide[group == "healthy"][[s]])
  rows[[length(rows) + 1]] <- data.table(comparison = paste("patient minus healthy,", s), n = nrow(wide),
                                         estimate = unname(t$estimate[1] - t$estimate[2]),
                                         lower = t$conf.int[1], upper = t$conf.int[2])
}
for (s in c("neb", "all10")) {
  shift <- wide[[s]] - wide[["none"]]
  t <- t.test(shift[wide$group == "patient"], shift[wide$group == "healthy"])
  rows[[length(rows) + 1]] <- data.table(comparison = paste("patient minus healthy, of", s, "minus none"), n = nrow(wide),
                                         estimate = unname(t$estimate[1] - t$estimate[2]),
                                         lower = t$conf.int[1], upper = t$conf.int[2])
}
fwrite(rbindlist(rows), file.path(out, "global_diff.tsv"), sep = "\t")

# ---- tiles ----
read_diff <- function(s) {
  f <- Sys.glob(file.path(work, "dmr", s, "methylDiff_*.txt.bgz"))
  stopifnot(length(f) == 1)
  # methylKit writes its object metadata as leading '#' lines in the tabix file
  d <- fread(cmd = paste("zcat", shQuote(f), "| grep -v '^#'"), header = FALSE,
             col.names = c("chr", "start", "end", "strand", "pvalue", "qvalue", "meth.diff"))
  d[, tile := paste(chr, start)]
  setkey(d, tile)
  d
}
diffs <- setNames(lapply(settings, read_diff), settings)
called <- function(d, q, delta) d[qvalue < q & abs(meth.diff) >= delta]
thresholds <- list(c(0.01, 25), c(0.05, 10))  # methylKit's getMethylDiff default, and a looser one
tiles <- rbindlist(lapply(settings, function(s) {
  d <- diffs[[s]]
  data.table(setting = s, tiles_tested = nrow(d),
             called_q01_d25 = nrow(called(d, 0.01, 25)), called_q05_d10 = nrow(called(d, 0.05, 10)))
}))
fwrite(tiles, file.path(out, "tiles.tsv"), sep = "\t")

pair_rows <- list()
for (p in pairs) {
  a <- diffs[[p[1]]]; b <- diffs[[p[2]]]
  m <- merge(a, b, by = "tile", suffixes = c(".a", ".b"))
  for (th in thresholds) {
    ca <- called(a, th[1], th[2])$tile; cb <- called(b, th[1], th[2])$tile
    both <- intersect(ca, cb)
    flips <- m[tile %in% both & sign(meth.diff.a) != sign(meth.diff.b), .N]
    pair_rows[[length(pair_rows) + 1]] <- data.table(
      comparison = paste(p[1], "against", p[2]), q = th[1], min_diff = th[2],
      tiles_shared = nrow(m), tiles_only_first = nrow(a) - nrow(m), tiles_only_second = nrow(b) - nrow(m),
      r_meth_diff = cor(m$meth.diff.a, m$meth.diff.b),
      median_abs_change = median(abs(m$meth.diff.a - m$meth.diff.b)),
      max_abs_change = max(abs(m$meth.diff.a - m$meth.diff.b)),
      called_first = length(ca), called_second = length(cb), called_both = length(both),
      jaccard = length(both) / length(union(ca, cb)), sign_differs = flips)
  }
}
fwrite(rbindlist(pair_rows), file.path(out, "tiles_pair.tsv"), sep = "\t")
