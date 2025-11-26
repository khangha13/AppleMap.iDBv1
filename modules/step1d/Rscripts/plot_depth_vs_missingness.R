#!/usr/bin/env Rscript

args_full <- commandArgs(trailingOnly = FALSE)
script_path <- sub(
  "^--file=", "",
  args_full[grepl("^--file=", args_full)]
)
script_dir <- if (length(script_path) > 0) {
  dirname(normalizePath(script_path[1]))
} else {
  normalizePath(".")
}
source(file.path(script_dir, "common_plot_utils.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop(
    "Usage: Rscript plot_depth_vs_missingness.R <metrics_tsv> <output_dir> [image_format]",
    call. = FALSE
  )
}

metrics_path <- normalizePath(args[1], mustWork = FALSE)
output_dir <- args[2]
image_format <- if (length(args) >= 3) args[3] else "png"

required_cols <- c("CHROM", "MEAN_DEPTH", "MISSING_RATE")
metrics_dt <- read_metrics_dataset(metrics_path, required_cols)

ensure_directory(output_dir)

chromosomes <- sort(unique(metrics_dt$CHROM))

for (chrom in chromosomes) {
  chr_data <- metrics_dt[
    CHROM == chrom &
      !is.na(MEAN_DEPTH) &
      !is.na(MISSING_RATE)
  ]
  if (nrow(chr_data) == 0) {
    next
  }

  plot_obj <- ggplot(chr_data, aes(x = MEAN_DEPTH, y = MISSING_RATE)) +
    geom_point(alpha = 0.35, color = "#6a3d9a", size = 0.5) +
    scale_x_continuous(limits = c(0, 150)) +
    scale_y_continuous(
      labels = percent_format(accuracy = 1)
    ) +
    labs(
      title = sprintf("Depth vs Missingness — %s", chrom),
      x = "Mean Depth",
      y = "Missingness"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )

  output_path <- file.path(
    output_dir,
    sprintf("%s_depth_vs_missingness.%s", chrom, image_format)
  )
  save_plot(
    plot_obj, output_path,
    width = 6, height = 6, device = image_format
  )
}