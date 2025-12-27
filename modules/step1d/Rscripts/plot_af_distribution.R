#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package data.table is required.", call. = FALSE)
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package ggplot2 is required.", call. = FALSE)
  }
  if (!requireNamespace("ragg", quietly = TRUE)) {
    stop("Package ragg is required.", call. = FALSE)
  }
  if (!requireNamespace("scales", quietly = TRUE)) {
    stop("Package scales is required.", call. = FALSE)
  }
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript plot_af_distribution.R <site_metrics.tsv> <output_dir> <image_format> [bins]", call. = FALSE)
}

metrics_path <- args[[1]]
output_dir <- args[[2]]
img_format <- args[[3]]
bins <- if (length(args) >= 4) as.integer(args[[4]]) else 50L
if (is.na(bins) || bins <= 0) bins <- 50L

if (!file.exists(metrics_path)) {
  stop(sprintf("Site metrics file not found: %s", metrics_path), call. = FALSE)
}
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

dt <- data.table::fread(metrics_path, data.table = FALSE)
if (!all(c("CHROM", "AF") %in% colnames(dt))) {
  stop("Site metrics file must contain CHROM and AF columns.", call. = FALSE)
}
dt$AF <- suppressWarnings(as.numeric(dt$AF))
dt <- dt[!is.na(dt$AF) & dt$AF >= 0 & dt$AF <= 1, , drop = FALSE]
if (nrow(dt) == 0) {
  stop("No valid allele frequency values found in the metrics file.", call. = FALSE)
}

library(ggplot2)

save_plot <- function(plot, path, fmt) {
  fmt <- tolower(fmt)
  if (fmt == "png") {
    ragg::agg_png(filename = path, width = 1600, height = 1200, res = 200)
    print(plot)
    dev.off()
  } else {
    ggplot2::ggsave(filename = path, plot = plot, width = 8, height = 6, device = fmt)
  }
}

make_hist <- function(df, title_suffix) {
  ggplot(df, aes(x = AF)) +
    geom_histogram(bins = bins, fill = "#1f78b4", colour = "white", alpha = 0.9) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1), labels = scales::number_format(accuracy = 0.1)) +
    theme_minimal(base_size = 14) +
    labs(
      title = sprintf("Allele Frequency Distribution%s", title_suffix),
      x = "Allele frequency",
      y = "Variant count"
    )
}

chroms <- sort(unique(dt$CHROM))
for (chr in chroms) {
  chr_df <- dt[dt$CHROM == chr, , drop = FALSE]
  if (nrow(chr_df) == 0) next
  plot <- make_hist(chr_df, sprintf(" – %s", chr))
  out_path <- file.path(output_dir, sprintf("%s_af_distribution.%s", chr, img_format))
  save_plot(plot, out_path, img_format)
}

# Combined plot (all chromosomes) for quick inspection
combined_plot <- make_hist(dt, "")
combined_path <- file.path(output_dir, sprintf("all_chromosomes_af_distribution.%s", img_format))
save_plot(combined_plot, combined_path, img_format)

message(sprintf("Allele frequency distribution plots written to %s", normalizePath(output_dir)))
