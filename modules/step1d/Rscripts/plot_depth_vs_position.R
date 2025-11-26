#!/usr/bin/env Rscript

args_full <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args_full[grepl("^--file=", args_full)])
script_dir <- if (length(script_path) > 0) dirname(normalizePath(script_path[1])) else normalizePath(".")
source(file.path(script_dir, "common_plot_utils.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript plot_depth_vs_position.R <metrics_tsv> <output_dir> [image_format]", call. = FALSE)
}

metrics_path <- normalizePath(args[1], mustWork = FALSE)
output_dir <- args[2]
image_format <- if (length(args) >= 3) args[3] else "png"

required_cols <- c("CHROM", "POS", "MEAN_DEPTH")
metrics_dt <- read_metrics_dataset(metrics_path, required_cols)
metrics_dt <- normalise_positions(metrics_dt)

ensure_directory(output_dir)

chromosomes <- sort(unique(metrics_dt$CHROM))

for (chrom in chromosomes) {
  chr_data <- metrics_dt[CHROM == chrom & !is.na(MEAN_DEPTH)]
  if (nrow(chr_data) == 0) {
    next
  }

  min_pos <- floor(min(chr_data$POS_MBP, na.rm = TRUE) / 5) * 5
  max_pos <- ceiling(max(chr_data$POS_MBP, na.rm = TRUE) / 5) * 5
  if (!is.finite(min_pos) || !is.finite(max_pos)) {
    next
  }
  if (max_pos <= min_pos) {
    max_pos <- min_pos + 5
  }
  breaks <- seq(min_pos, max_pos, by = 5)

  plot_obj <- ggplot(chr_data, aes(x = POS_MBP, y = MEAN_DEPTH)) +
    geom_line(color = "#1f78b4", linewidth = 0.5) +
    geom_point(color = "#1f78b4", alpha = 0.4, size = 0.4) +
    scale_x_continuous(limits = c(min_pos, max_pos), breaks = breaks) +
    labs(
      title = sprintf("Mean Depth Across %s", chrom),
      x = "Position (Mbp)",
      y = "Mean Depth"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
  
  output_path <- file.path(output_dir, sprintf("%s_depth_vs_position.%s", chrom, image_format))
  save_plot(plot_obj, output_path, width = 12, height = 4, device = image_format)
}
