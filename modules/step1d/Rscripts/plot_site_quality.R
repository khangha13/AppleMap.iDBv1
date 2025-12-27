#!/usr/bin/env Rscript

args_full <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args_full[grepl("^--file=", args_full)])
script_dir <- if (length(script_path) > 0) dirname(normalizePath(script_path[1])) else normalizePath(".")
source(file.path(script_dir, "common_plot_utils.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript plot_site_quality.R <metrics_tsv> <output_dir> [image_format]", call. = FALSE)
}

metrics_path <- normalizePath(args[1], mustWork = FALSE)
output_dir <- args[2]
image_format <- if (length(args) >= 3) args[3] else "png"

required_cols <- c("CHROM", "POS", "QUAL")
metrics_dt <- read_metrics_dataset(metrics_path, required_cols)
metrics_dt <- normalise_positions(metrics_dt)

ensure_directory(output_dir)

chromosomes <- sort(unique(metrics_dt$CHROM))

for (chrom in chromosomes) {
  chr_data <- metrics_dt[CHROM == chrom & !is.na(QUAL)]
  if (nrow(chr_data) == 0) {
    next
  }

  # Ensure QUAL is numeric (Phred scores from VCF)
  chr_data[, QUAL := as.numeric(QUAL)]
  
  # Filter out extreme outliers (QUAL > 100 is unrealistic for Phred scores)
  # Cap at 100 to focus on meaningful quality range
  chr_data <- chr_data[is.finite(QUAL) & QUAL >= 0 & QUAL <= 100]
  
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

  plot_obj <- ggplot(chr_data, aes(x = POS_MBP, y = QUAL)) +
    geom_point(alpha = 0.35, color = "#ff7f00", size = 0.5) +
    scale_x_continuous(limits = c(min_pos, max_pos), breaks = breaks) +
    scale_y_continuous(limits = c(0, 100)) +
    labs(
      title = sprintf("Site Quality Across %s", chrom),
      x = "Position (Mbp)",
      y = "QUAL (Phred Score)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
  
  output_path <- file.path(output_dir, sprintf("%s_site_quality.%s", chrom, image_format))
  save_plot(plot_obj, output_path, width = 12, height = 4, device = image_format)
}
