#!/usr/bin/env Rscript

args_full <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args_full[grepl("^--file=", args_full)])
script_dir <- if (length(script_path) > 0) dirname(normalizePath(script_path[1])) else normalizePath(".")
source(file.path(script_dir, "common_plot_utils.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript plot_quality_by_depth.R <metrics_tsv> <output_dir> [image_format] [bins]", call. = FALSE)
}

metrics_path <- normalizePath(args[1], mustWork = FALSE)
output_dir <- args[2]
image_format <- if (length(args) >= 3) args[3] else "png"
bins <- if (length(args) >= 4) as.integer(args[4]) else 40L
if (is.na(bins) || bins <= 0) {
  stop("Bins parameter must be a positive integer.", call. = FALSE)
}

required_cols <- c("CHROM", "QD")
metrics_dt <- read_metrics_dataset(metrics_path, required_cols)

metrics_dt[, QD := as.numeric(QD)]
metrics_dt <- metrics_dt[!is.na(QD) & is.finite(QD)]

if (nrow(metrics_dt) == 0) {
  stop("No valid QD values found in metrics dataset.", call. = FALSE)
}

metrics_dt[, CHROM := factor(CHROM, levels = unique(metrics_dt$CHROM))]

ensure_directory(output_dir)

plot_obj <- ggplot(metrics_dt, aes(x = QD)) +
  geom_histogram(bins = bins, fill = "#1c9099", color = "#0c4a63", alpha = 0.85) +
  facet_wrap(~ CHROM, scales = "free_y") +
  labs(
    title = "Quality-by-Depth (QD) Distribution",
    x = "QD",
    y = "Variant Count"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold")
  )

output_path <- file.path(output_dir, sprintf("quality_by_depth_hist.%s", image_format))
save_plot(plot_obj, output_path, width = 12, height = 7, device = image_format)
