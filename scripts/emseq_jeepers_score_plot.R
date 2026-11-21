#!/usr/bin/env Rscript
# Figure of the JEEPERS test on one plasma library.
#
# Usage: emseq_jeepers_score_plot.R <results_dir> <out.png> [figure_schema.R]
#
# Reads gap.tsv and opposite_strand.tsv written by emseq_jeepers_score.py.
# Left: CpG methylation by distance from the 3' end and from the 5' end of the
# strand, on fragments of 160 bp or longer. Right: the share of a strand's
# calls that equal the call of the opposite strand of the same molecule, by
# distance from the strand's 3' end. Points are shares with 95% Wilson
# intervals. Uncorrected calls and the locus majority call are drawn from the
# rows of the code-default setting. The image carries no caption; the report
# supplies it.

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: emseq_jeepers_score_plot.R <results_dir> <out.png> [figure_schema.R]")
schema <- if (length(args) >= 3) args[3] else "~/repos/science/R/figure_schema.R"
source(schema)

bins <- c("1-10", "11-20", "21-30", "31-50", "51-80")
series <- c(
  uncorrected = "Uncorrected",
  jeepers_0.9 = "JEEPERS, -m 0.9",
  jeepers_default = "JEEPERS, code default",
  majority = "Locus majority call",
  excl30 = "Exclude 30 bases",
  excl52 = "Exclude 52 bases"
)
palette <- setNames(c("#4D4D4D", "#D55E00", "#E69F00", "#0072B2", "#009E73", "#CC79A7"), series)

# one series label per row: JEEPERS from both settings, the rest from the code-default setting
label_rows <- function(d) {
  d$key <- ifelse(d$method == "jeepers", paste0("jeepers_", d$setting), d$method)
  d <- d[d$method == "jeepers" | d$setting == "default", ]
  d <- d[d$key %in% names(series), ]
  d$series <- factor(series[d$key], levels = series)
  d
}

wilson <- function(k, n) {
  z <- qnorm(0.975)
  p <- k / n
  mid <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
  data.frame(pct = 100 * p, lower = 100 * (mid - half), upper = 100 * (mid + half))
}

gap <- read.delim(file.path(args[1], "gap.tsv"), stringsAsFactors = FALSE, colClasses = c(setting = "character"))
gap <- label_rows(gap)
ends <- c("3p" = "3' end of strand", "5p" = "5' end of strand")
meth <- do.call(rbind, lapply(names(ends), function(e) {
  n <- gap[[paste0("n_", e)]]
  cbind(gap[, c("series", "bin")], end = ends[[e]], n = n, wilson(gap[[paste0("meth_", e)]], n))
}))
meth <- meth[meth$n > 0, ]
meth$end <- factor(meth$end, levels = ends)
meth$bin <- factor(meth$bin, levels = bins)

opp <- read.delim(file.path(args[1], "opposite_strand.tsv"), stringsAsFactors = FALSE,
                  colClasses = c(setting = "character"))
opp <- opp[opp$pairs == "same" & opp$row == "agreement" & opp$scope %in% bins, ]
opp <- label_rows(opp)
agree <- cbind(opp[, c("series", "scope", "n")], wilson(opp$count, opp$n))
agree <- agree[agree$n > 0, ]
agree$bin <- factor(agree$scope, levels = bins)

dodge <- position_dodge(width = 0.6)
layers <- list(
  geom_line(aes(group = series), position = dodge, linewidth = 0.4),
  geom_errorbar(aes(ymin = lower, ymax = upper), position = dodge, width = 0, linewidth = 0.4),
  geom_point(position = dodge, size = 1.4),
  scale_color_manual(values = palette, name = NULL, drop = FALSE),
  scale_x_discrete(drop = FALSE),
  theme_scifig(base_size = 11),
  theme(axis.text = element_text(color = "black"), strip.text = element_text(hjust = 0))
)

p_meth <- ggplot(meth, aes(bin, pct, color = series)) + layers +
  facet_wrap(~end, nrow = 1) +
  labs(x = "Distance from strand end (bases)", y = "CpG methylation (%)", tag = "A") +
  theme(legend.position = "none")

p_agree <- ggplot(agree, aes(bin, pct, color = series)) + layers +
  labs(x = "Distance from 3' end of strand (bases)", y = "Agreement with opposite strand (%)", tag = "B") +
  theme(legend.position = "right", legend.key.height = unit(0.5, "cm"))

png(args[2], width = 12.5, height = 3.8, units = "in", res = 300, bg = "white")
pushViewport(viewport(layout = grid.layout(1, 2, widths = unit(c(1.1, 1), "null"))))
print(p_meth, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p_agree, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
invisible(dev.off())
cat("Saved:", args[2], "\n")
