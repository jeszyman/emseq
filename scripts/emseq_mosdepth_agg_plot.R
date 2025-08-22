#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  suppressWarnings(library(argparse))
  suppressWarnings(library(data.table))
  suppressWarnings(library(ggplot2))
  suppressWarnings(library(Cairo))
  suppressWarnings(library(scales))
  suppressWarnings(library(patchwork))
  suppressWarnings(library(matrixStats))
})

# -------------------------------
# Argument parsing
# -------------------------------

prog <- basename(commandArgs(trailingOnly = FALSE)[1])

parser <- ArgumentParser(
  description = "Generate a paginated threshold coverage plot from mosdepth output.",
  prog = prog
)

parser$add_argument("--threshold_list", required = TRUE,
                    help = "Space-separated list of mosdepth threshold files (*.thresholds.bed.gz)")
parser$add_argument("--regions_list", required = TRUE,
                    help = "Space-separated list of mosdepth regions files (*.regions.bed.gz)")
parser$add_argument("--library_list", required = TRUE,
                    help = "Space-separated list of sample names (must match file order)")
parser$add_argument("--output_pdf", required = TRUE,
                    help = "Full output PDF file path (e.g., /tmp/plot.pdf)")
parser$add_argument("--output_tsv", required = TRUE,
                    help = "Path for tabular output")

args <- parser$parse_args()
threshold_files <- unlist(strsplit(args$threshold_list, " "))
regions_files <- unlist(strsplit(args$regions_list, " "))
library_ids <- unlist(strsplit(args$library_list, " "))
output_tsv <- args$output_tsv
output_pdf <- args$output_pdf

if (!all(lengths(list(threshold_files, regions_files, library_ids)) == length(library_ids))) {
  stop("Error: threshold_list, regions_list, and library_list must all be the same length.")
}

# -------------------------------
# Read and melt threshold files
# -------------------------------

read_thresholds <- function(file, sample) {
  header <- fread(file, nrows = 0)
  names(header)[1] <- sub("^#", "", names(header)[1])
  threshold_cols <- setdiff(names(header), c("chrom", "start", "end", "region"))
  df <- fread(file, skip = 1, col.names = names(header))
  df[, sample := sample]
  melted <- melt(df,
    id.vars = c("chrom", "start", "end", "region", "sample"),
    measure.vars = threshold_cols,
    variable.name = "threshold",
    value.name = "count"
  )
  list(data = melted, thresholds = threshold_cols)
}

parsed <- mapply(read_thresholds, threshold_files, library_ids, SIMPLIFY = FALSE)
hist_data <- rbindlist(lapply(parsed, `[[`, "data"))
hist_data[, count := as.numeric(count)]

all_thresholds <- unique(unlist(lapply(parsed, `[[`, "thresholds")))

# -------------------------------
# Read autosomal median from regions.bed.gz
# -------------------------------

get_autosomal_median <- function(file, sample_id) {
  df <- fread(file, col.names = c("chrom", "start", "end", "depth"))
  df <- df[chrom %in% paste0("chr", 1:22)]
  df[, sample := sample_id]
  df[, median := median(depth)]
  df[1, .(sample, median)]
}

medians <- rbindlist(mapply(get_autosomal_median, regions_files, library_ids, SIMPLIFY = FALSE))

# -------------------------------
# Infer 0X bins from zeroed rows
# -------------------------------

hist_wide <- dcast(hist_data, chrom + start + end + region + sample ~ threshold,
                   value.var = "count", fill = 0)
hist_wide[, is_zero := rowSums(.SD) == 0, .SDcols = all_thresholds]
zero_counts <- hist_data[, .(total = sum(count)), by = .(chrom, start, end, region, sample)]
zero_counts <- zero_counts[total == 0, .(count = .N * (end[1] - start[1])), by = sample]
zero_counts[, threshold := "0X"]

# -------------------------------
# Aggregate and bind all data
# -------------------------------

plot_data <- hist_data[, .(count = sum(count)), by = .(sample, threshold)]
plot_data <- rbind(plot_data, zero_counts, fill = TRUE)

threshold_levels <- unique(plot_data$threshold)
threshold_levels <- threshold_levels[order(as.numeric(sub("X$", "", as.character(threshold_levels))))]
plot_data[, threshold := factor(threshold, levels = threshold_levels)]

# -------------------------------
# Plot panels
# -------------------------------
print(medians)

make_panel <- function(sample_id) {
  median_val <- medians[sample == sample_id][["median"]]
  subtitle <- sprintf("Median depth: %.1f×", median_val)

  ggplot(plot_data[sample == sample_id], aes(x = threshold, y = count, fill = threshold)) +
    geom_col(width = 0.8) +
    scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
    scale_fill_brewer(palette = "Set2", guide = "none") +
    labs(title = sample_id, subtitle = subtitle, x = "Coverage threshold", y = "Covered bases") +
    theme_minimal(base_size = 10) +
    theme(
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 9),
      plot.title = element_text(size = 10, hjust = 0.5),
      plot.subtitle = element_text(size = 9, hjust = 0.5),
      panel.grid = element_line(linewidth = 0.2, colour = "grey90")
    )
}

ncol <- 4
nrow <- 6
panels_per_page <- ncol * nrow
sample_list <- unique(plot_data$sample)
pages <- split(sample_list, ceiling(seq_along(sample_list) / panels_per_page))

# -------------------------------
# Output
# -------------------------------

CairoPDF(output_pdf, width = 8.5, height = 11, onefile = TRUE)
for (i in seq_along(pages)) {
  plots <- lapply(pages[[i]], make_panel)
  layout <- wrap_plots(plots, ncol = ncol, nrow = nrow) +
    plot_annotation(
      title = "Coverage threshold by sample",
      theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
    )
  print(layout)
}
dev.off()

fwrite(medians,   file = args$output_tsv, sep = "\t")

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  suppressWarnings(library(argparse))
  suppressWarnings(library(data.table))
  suppressWarnings(library(ggplot2))
  suppressWarnings(library(Cairo))
  suppressWarnings(library(scales))
  suppressWarnings(library(patchwork))
  suppressWarnings(library(matrixStats))
})

# -------------------------------
# Argument parsing
# -------------------------------

prog <- basename(commandArgs(trailingOnly = FALSE)[1])

parser <- ArgumentParser(
  description = "Generate a paginated threshold coverage plot from mosdepth output.",
  prog = prog
)

parser$add_argument("--threshold_list", required = TRUE,
                    help = "Space-separated list of mosdepth threshold files (*.thresholds.bed.gz)")
parser$add_argument("--regions_list", required = TRUE,
                    help = "Space-separated list of mosdepth regions files (*.regions.bed.gz)")
parser$add_argument("--library_list", required = TRUE,
                    help = "Space-separated list of sample names (must match file order)")
parser$add_argument("--output_pdf", required = TRUE,
                    help = "Full output PDF file path (e.g., /tmp/plot.pdf)")
parser$add_argument("--output_tsv", required = TRUE,
                    help = "Path for tabular output")

args <- parser$parse_args()
threshold_files <- unlist(strsplit(args$threshold_list, " "))
regions_files <- unlist(strsplit(args$regions_list, " "))
library_ids <- unlist(strsplit(args$library_list, " "))
output_tsv <- args$output_tsv
output_pdf <- args$output_pdf

if (!all(lengths(list(threshold_files, regions_files, library_ids)) == length(library_ids))) {
  stop("Error: threshold_list, regions_list, and library_list must all be the same length.")
}

# -------------------------------
# Read and melt threshold files
# -------------------------------

read_thresholds <- function(file, sample) {
  header <- fread(file, nrows = 0)
  names(header)[1] <- sub("^#", "", names(header)[1])
  threshold_cols <- setdiff(names(header), c("chrom", "start", "end", "region"))
  df <- fread(file, skip = 1, col.names = names(header))
  df[, sample := sample]
  melted <- melt(df,
    id.vars = c("chrom", "start", "end", "region", "sample"),
    measure.vars = threshold_cols,
    variable.name = "threshold",
    value.name = "count"
  )
  list(data = melted, thresholds = threshold_cols)
}

parsed <- mapply(read_thresholds, threshold_files, library_ids, SIMPLIFY = FALSE)
hist_data <- rbindlist(lapply(parsed, `[[`, "data"))
hist_data[, count := as.numeric(count)]

all_thresholds <- unique(unlist(lapply(parsed, `[[`, "thresholds")))

# -------------------------------
# Read autosomal median from regions.bed.gz
# -------------------------------

get_autosomal_median <- function(file, sample_id) {
  df <- fread(file, col.names = c("chrom", "start", "end", "depth"))
  df <- df[chrom %in% paste0("chr", 1:22)]
  df[, sample := sample_id]
  df[, median := median(depth)]
  df[1, .(sample, median)]
}

medians <- rbindlist(mapply(get_autosomal_median, regions_files, library_ids, SIMPLIFY = FALSE))

# -------------------------------
# Infer 0X bins from zeroed rows
# -------------------------------

hist_wide <- dcast(hist_data, chrom + start + end + region + sample ~ threshold,
                   value.var = "count", fill = 0)
hist_wide[, is_zero := rowSums(.SD) == 0, .SDcols = all_thresholds]
zero_counts <- hist_data[, .(total = sum(count)), by = .(chrom, start, end, region, sample)]
zero_counts <- zero_counts[total == 0, .(count = .N * (end[1] - start[1])), by = sample]
zero_counts[, threshold := "0X"]

# -------------------------------
# Aggregate and bind all data
# -------------------------------

plot_data <- hist_data[, .(count = sum(count)), by = .(sample, threshold)]
plot_data <- rbind(plot_data, zero_counts, fill = TRUE)

threshold_levels <- unique(plot_data$threshold)
threshold_levels <- threshold_levels[order(as.numeric(sub("X$", "", as.character(threshold_levels))))]
plot_data[, threshold := factor(threshold, levels = threshold_levels)]

# -------------------------------
# Plot panels
# -------------------------------
print(medians)

make_panel <- function(sample_id) {
  median_val <- medians[sample == sample_id][["median"]]
  subtitle <- sprintf("Median depth: %.1f×", median_val)

  ggplot(plot_data[sample == sample_id], aes(x = threshold, y = count, fill = threshold)) +
    geom_col(width = 0.8) +
    scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
    scale_fill_brewer(palette = "Set2", guide = "none") +
    labs(title = sample_id, subtitle = subtitle, x = "Coverage threshold", y = "Covered bases") +
    theme_minimal(base_size = 10) +
    theme(
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 9),
      plot.title = element_text(size = 10, hjust = 0.5),
      plot.subtitle = element_text(size = 9, hjust = 0.5),
      panel.grid = element_line(linewidth = 0.2, colour = "grey90")
    )
}

ncol <- 4
nrow <- 6
panels_per_page <- ncol * nrow
sample_list <- unique(plot_data$sample)
pages <- split(sample_list, ceiling(seq_along(sample_list) / panels_per_page))

# -------------------------------
# Output
# -------------------------------

CairoPDF(output_pdf, width = 8.5, height = 11, onefile = TRUE)
for (i in seq_along(pages)) {
  plots <- lapply(pages[[i]], make_panel)
  layout <- wrap_plots(plots, ncol = ncol, nrow = nrow) +
    plot_annotation(
      title = "Coverage threshold by sample",
      theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
    )
  print(layout)
}
dev.off()

fwrite(medians,   file = args$output_tsv, sep = "\t")
