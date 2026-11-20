#!/usr/bin/env Rscript
# Control methylation by position from each read end (M-bias plot).
#
# Reads the position table of the read-end validation and draws CpG
# methylation of the pUC19 and lambda controls over the first bases from
# each of the four read ends, pooled over libraries, with 95% Wilson
# intervals. The image carries no caption; the report supplies it.

suppressPackageStartupMessages({
  library(argparse)
  library(ggplot2)
})

parser <- ArgumentParser(description = "M-bias plot of the lambda and pUC19 controls")
parser$add_argument("--position", required = TRUE, help = "position.tsv from the values script")
parser$add_argument("--values", required = TRUE, help = "values.json from the values script")
parser$add_argument("--out", required = TRUE, help = "output PNG")
parser$add_argument("--max-position", type = "integer", default = 30,
                    help = "bases from the read end to show")
args <- parser$parse_args()

# Strand whose table position 1 sits at each read end; MethylDackel counts
# positions from the left end of the stored alignment
# (tests/unit/test_methyldackel_bounds.py).
#
# Read 1 starts at the 5' end of the original strand and read 2 starts at its
# 3' end. The last base of a read is the opposite strand end only on a
# fragment the read runs through, so the two "end of read" panels mix strand
# ends with interior bases. Panels are ordered 5' side first.
ends <- data.frame(
  end    = c("5' end of strand\nStart of read 1", "5' end of strand\nEnd of read 2",
             "3' end of strand\nStart of read 2", "3' end of strand\nEnd of read 1"),
  strand = c("OT", "OT", "OB", "OB"),
  read   = c(1, 2, 2, 1)
)
controls <- c(puc19 = "pUC19 (CpG methylated)", unmeth_lambda = "Lambda (unmethylated)")

d <- read.delim(args$position, stringsAsFactors = FALSE)
d <- d[d$library_index == 0 & d$context == "CpG" & d$position <= args$max_position, ]
d <- merge(d, ends, by = c("strand", "read"))
d$end <- factor(d$end, levels = ends$end)
d$control <- factor(controls[d$source], levels = controls)
d$pct <- 100 * d$n_methylated / (d$n_methylated + d$n_unmethylated)
d$lower_pct <- 100 * d$lower
d$upper_pct <- 100 * d$upper

palette <- setNames(c("#1F5F99", "#D9730D"), controls)

p <- ggplot(d, aes(position, pct, color = control, fill = control)) +
  geom_ribbon(aes(ymin = lower_pct, ymax = upper_pct), alpha = 0.25, color = NA) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1.1) +
  facet_wrap(~end, nrow = 1) +
  scale_color_manual(values = palette, name = NULL) +
  scale_fill_manual(values = palette, name = NULL) +
  scale_x_continuous(breaks = c(1, 10, 20, 30), expand = expansion(mult = 0.03)) +
  scale_y_continuous(breaks = seq(0, 100, 25), limits = c(0, 100), expand = expansion(mult = 0.03)) +
  labs(x = "Bases from read end", y = "CpG methylation (%)") +
  theme_classic(base_size = 11, base_family = "sans") +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 11, hjust = 0),
    axis.title.x = element_text(margin = margin(t = 8)),
    axis.title.y = element_text(margin = margin(r = 8)),
    axis.text = element_text(color = "black"),
    panel.spacing.x = unit(14, "pt"),
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.4, "cm"),
    legend.margin = margin(t = 0),
    plot.margin = margin(8, 12, 6, 8, "pt")
  )

ggsave(args$out, p, width = 8.5, height = 3.4, dpi = 300, bg = "white")
