#!/usr/bin/env Rscript

# Shared helper functions for Step 1D plotting scripts

suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required but not installed.", call. = FALSE)
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required but not installed.", call. = FALSE)
  }
  if (!requireNamespace("scales", quietly = TRUE)) {
    stop("Package 'scales' is required but not installed.", call. = FALSE)
  }
})

library(data.table)
library(ggplot2)
library(scales)

ensure_directory <- function(dir_path) {
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  }
}

read_metrics_dataset <- function(metrics_path, required_cols) {
  if (!file.exists(metrics_path)) {
    stop(sprintf("Metrics file not found: %s", metrics_path), call. = FALSE)
  }
  dt <- data.table::fread(metrics_path, na.strings = c("NA", "NaN", ""))
  missing_cols <- setdiff(required_cols, names(dt))
  if (length(missing_cols) > 0) {
    stop(
      sprintf(
        "Metrics file is missing required columns: %s",
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  numeric_candidates <- c(
    "POS", "QUAL", "QD", "AC", "AF", "INBREEDING_COEFF", "EXCESS_HET", "MQ",
    "MEAN_DEPTH", "CALL_RATE", "MISSING_RATE", "HETEROZYGOUS_RATE", "DR2", "IMP",
    "DP_NON_MISSING", "CALLED_GENOTYPES", "MISSING_GENOTYPES",
    "TOTAL_GENOTYPES", "HETEROZYGOUS_COUNT"
  )
  numeric_cols <- intersect(numeric_candidates, names(dt))
  for (col in numeric_cols) {
    dt[[col]] <- as.numeric(dt[[col]])
  }
  unique(dt, by = c("CHROM", "POS"))
}

normalise_positions <- function(dt) {
  dt[, POS_MBP := POS / 1e6]
  dt
}

save_plot <- function(plot_obj, output_path, width = 12, height = 4, dpi = 300, device = NULL) {
  if (is.null(device)) {
    device <- tools::file_ext(output_path)
  }
  ggplot2::ggsave(
    filename = output_path,
    plot = plot_obj,
    width = width,
    height = height,
    dpi = dpi,
    units = "in",
    device = device,
    limitsize = FALSE
  )
}
