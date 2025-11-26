#!/usr/bin/env Rscript

args_full <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args_full[grepl("^--file=", args_full)])
script_dir <- if (length(script_path) > 0) dirname(normalizePath(script_path[1])) else normalizePath(".")
source(file.path(script_dir, "common_plot_utils.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript plot_call_rate_heatmap.R <metrics_tsv> <output_dir> [image_format] [bins]", call. = FALSE)
}

metrics_path <- normalizePath(args[1], mustWork = FALSE)
output_dir <- args[2]
image_format <- if (length(args) >= 3) args[3] else "png"
bins <- if (length(args) >= 4) as.integer(args[4]) else 100L
if (is.na(bins) || bins <= 0) {
  stop("Bins parameter must be a positive integer.", call. = FALSE)
}

required_cols <- c("CHROM", "POS", "CALL_RATE")
metrics_dt <- read_metrics_dataset(metrics_path, required_cols)
metrics_dt <- normalise_positions(metrics_dt)

ensure_directory(output_dir)

chromosomes <- sort(unique(metrics_dt$CHROM))

for (chrom in chromosomes) {
  chr_data <- metrics_dt[CHROM == chrom & !is.na(CALL_RATE)]
  if (nrow(chr_data) == 0) {
    next
  }
  
  if (nrow(chr_data) < 5) {
    warning(sprintf("Skipping %s heatmap — insufficient variant sites.", chrom))
    next
  }
  
  min_pos <- min(chr_data$POS_MBP, na.rm = TRUE)
  max_pos <- max(chr_data$POS_MBP, na.rm = TRUE)
  
  if (!is.finite(min_pos) || !is.finite(max_pos)) {
    next
  }
  
  if (max_pos == min_pos) {
    agg <- data.table(
      BIN = 1L,
      CALL_RATE = mean(chr_data$CALL_RATE, na.rm = TRUE),
      BIN_CENTER = min_pos
    )
  } else {
    breaks <- seq(min_pos, max_pos, length.out = bins + 1)
    chr_data[, BIN := cut(POS_MBP, breaks = breaks, include.lowest = TRUE, labels = FALSE)]
    
    agg <- chr_data[!is.na(BIN), .(
      CALL_RATE = mean(CALL_RATE, na.rm = TRUE),
      BIN_CENTER = mean(POS_MBP, na.rm = TRUE)
    ), by = BIN]
  }
  
  if (nrow(agg) == 0) {
    warning(sprintf("Skipping %s heatmap — unable to compute bin aggregates.", chrom))
    next
  }
  
  plot_obj <- ggplot(agg, aes(x = BIN_CENTER, y = 1, fill = CALL_RATE)) +
    geom_tile(color = NA, height = 1) +
    scale_fill_gradientn(
      colours = c("#b2182b", "#fddbc7", "#d1e5f0", "#2166ac"),
      limits = c(0, 1),
      labels = percent_format(accuracy = 1),
      name = "Call Rate"
    ) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(
      title = sprintf("Call Rate Heat Map — %s", chrom),
      x = "Position (Mbp)",
      y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      legend.position = "bottom"
    )
  
  output_path <- file.path(output_dir, sprintf("%s_call_rate_heatmap.%s", chrom, image_format))
  save_plot(plot_obj, output_path, width = 12, height = 2.5, device = image_format)
}

