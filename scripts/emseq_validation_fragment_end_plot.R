#!/usr/bin/env Rscript
# CpG methylation by distance from the 5' and 3' ends of the original strand.
#
# Reads one or more fragment-end tables of the read-end validation and draws
# one panel per table: methylation over the first bases from each strand end,
# on fragments long enough that a base is nearer its own end than the other,
# counted once per fragment and pooled over the judged libraries, with 95%
# Wilson intervals. Single libraries are drawn as thin lines behind the pooled
# line. The image carries no caption; the report supplies it.

suppressPackageStartupMessages({
  library(argparse)
  library(ggplot2)
})

parser <- ArgumentParser(description = "Methylation by distance from the strand ends")
parser$add_argument("--table", required = TRUE, nargs = "+", help = "fragment_end.tsv files, one per panel")
parser$add_argument("--source", required = TRUE, nargs = "+", help = "source to draw from each table")
parser$add_argument("--label", required = TRUE, nargs = "+", help = "panel title for each table")
parser$add_argument("--exclude-index", type = "integer", nargs = "*",
                    help = "library index to leave out of each panel, 0 for none")
parser$add_argument("--out", required = TRUE, help = "output PNG")
args <- parser$parse_args()

end_labels <- c("5p" = "5' end of strand", "3p" = "3' end of strand")
z <- qnorm(0.975)

panels <- lapply(seq_along(args$table), function(i) {
  d <- read.delim(args$table[i], stringsAsFactors = FALSE)
  d <- d[d$source == args$source[i] & d$fragments == "long" & d$mate == 0 & d$distance <= 80, ]
  excluded <- if (!is.null(args$exclude_index) && length(args$exclude_index) >= i) args$exclude_index[i] else 0
  d <- d[d$library_index != excluded | excluded == 0, ]
  d$panel <- args$label[i]
  d
})
d <- do.call(rbind, panels)
d$panel <- factor(d$panel, levels = args$label)
d$end <- factor(end_labels[d$end], levels = end_labels)
d$n <- d$n_methylated + d$n_unmethylated
d <- d[d$n > 0, ]
prop <- d$n_methylated / d$n
mid <- (prop + z^2 / (2 * d$n)) / (1 + z^2 / d$n)
half <- z * sqrt(prop * (1 - prop) / d$n + z^2 / (4 * d$n^2)) / (1 + z^2 / d$n)
d$pct <- 100 * prop
d$lower_pct <- 100 * (mid - half)
d$upper_pct <- 100 * (mid + half)

pooled <- d[d$library_index == 0, ]
single <- d[d$library_index != 0, ]
palette <- setNames(c("#1F5F99", "#D9730D"), end_labels)

p <- ggplot(pooled, aes(distance, pct, color = end, fill = end)) +
  geom_line(data = single, aes(group = interaction(library_index, end)), linewidth = 0.2, alpha = 0.45) +
  geom_ribbon(aes(ymin = lower_pct, ymax = upper_pct), alpha = 0.25, color = NA) +
  geom_line(linewidth = 0.6) +
  facet_wrap(~panel, nrow = 1) +
  scale_color_manual(values = palette, name = NULL) +
  scale_fill_manual(values = palette, name = NULL) +
  scale_x_continuous(breaks = c(1, 20, 40, 60, 80), expand = expansion(mult = 0.02)) +
  scale_y_continuous(breaks = seq(0, 100, 25), limits = c(0, 100), expand = expansion(mult = 0.03)) +
  labs(x = "Bases from strand end", y = "CpG methylation (%)") +
  theme_classic(base_size = 11, base_family = "sans") +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 11, hjust = 0),
    axis.title.x = element_text(margin = margin(t = 8)),
    axis.title.y = element_text(margin = margin(r = 8)),
    axis.text = element_text(color = "black"),
    panel.spacing.x = unit(18, "pt"),
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.4, "cm"),
    legend.margin = margin(t = 0),
    plot.margin = margin(8, 12, 6, 8, "pt")
  )

ggsave(args$out, p, width = 3.4 * length(args$table) + 0.6, height = 3.6, dpi = 300, bg = "white")
