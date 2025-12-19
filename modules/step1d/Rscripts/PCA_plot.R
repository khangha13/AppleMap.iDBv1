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
  # Optional: ggrepel for nicer, non-overlapping text labels
  has_ggrepel <<- requireNamespace("ggrepel", quietly = TRUE)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript PCA_plot.R <eigenvec> <eigenval> <output_dir> [show_labels] [label_size] [use_ggrepel] [duplicate_samples]", call. = FALSE)
}

eigenvec_path <- args[[1]]
eigenval_path <- args[[2]]
output_dir <- args[[3]]
show_labels <- if (length(args) >= 4) as.logical(args[[4]]) else TRUE
label_size <- if (length(args) >= 5) as.numeric(args[[5]]) else 3
use_ggrepel <- if (length(args) >= 6) as.logical(args[[6]]) else TRUE
duplicate_samples_path <- if (length(args) >= 7) args[[7]] else ""

assert_file <- function(path, label) {
  if (!file.exists(path)) {
    stop(sprintf("%s not found: %s", label, path), call. = FALSE)
  }
}

assert_file(eigenvec_path, "Eigenvector file")
assert_file(eigenval_path, "Eigenvalue file")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

eig <- data.table::fread(eigenvec_path, data.table = FALSE)
if (ncol(eig) < 4) {
  stop("Eigenvec file must contain at least two PCs (PC1 & PC2).", call. = FALSE)
}
colnames(eig) <- c("FID", "IID", paste0("PC", seq_len(ncol(eig) - 2)))

dup_ids <- character(0)
if (nzchar(duplicate_samples_path) && file.exists(duplicate_samples_path)) {
  if (file.info(duplicate_samples_path)$size > 0) {
    dup_tbl <- data.table::fread(duplicate_samples_path, header = FALSE, data.table = FALSE)
    if (ncol(dup_tbl) == 1) {
      dup_ids <- dup_tbl[[1]]
    } else if (ncol(dup_tbl) >= 2) {
      dup_ids <- dup_tbl[[2]]
    }
    dup_ids <- unique(dup_ids)
  }
}

eigenvalues <- scan(eigenval_path)
if (length(eigenvalues) == 0) {
  stop("Eigenvalue file is empty.", call. = FALSE)
}
variance_prop <- eigenvalues / sum(eigenvalues)

axis_label <- function(pc_index) {
  pct <- if (pc_index <= length(variance_prop)) variance_prop[pc_index] * 100 else NA_real_
  if (is.na(pct)) {
    sprintf("PC%d", pc_index)
  } else {
    sprintf("PC%d (%.2f%%)", pc_index, pct)
  }
}

library(ggplot2)

dup_suffix <- ""
if (length(dup_ids) > 0) {
  dup_suffix <- sprintf(" • Duplicates flagged: %d", length(dup_ids))
}

if (length(dup_ids) > 0) {
  eig$DupStatus <- ifelse(eig$IID %in% dup_ids, "Duplicate", "Unique")
  scatter <- ggplot(eig, aes(x = PC1, y = PC2, colour = DupStatus)) +
    geom_point(alpha = 0.6, size = 0.9) +
    scale_color_manual(values = c("Unique" = "#1f78b4", "Duplicate" = "#e31a1c")) +
    theme_minimal(base_size = 14) +
    labs(
      title = "Principal Component Analysis",
      subtitle = sprintf("Samples: %d • Variants: %d (after pruning)%s", nrow(eig), length(eigenvalues), dup_suffix),
      x = axis_label(1),
      y = axis_label(2),
      colour = "Sample status"
    )
} else {
  scatter <- ggplot(eig, aes(x = PC1, y = PC2)) +
    geom_point(alpha = 0.6, size = 0.8, colour = "#1f78b4") +
    theme_minimal(base_size = 14) +
    labs(
      title = "Principal Component Analysis",
      subtitle = sprintf("Samples: %d • Variants: %d (after pruning)", nrow(eig), length(eigenvalues)),
      x = axis_label(1),
      y = axis_label(2)
    )
}

if (isTRUE(show_labels)) {
  if (isTRUE(use_ggrepel) && isTRUE(has_ggrepel)) {
    # Use ggrepel to avoid overlapping labels
    scatter <- scatter + ggrepel::geom_text_repel(
      aes(label = IID),
      colour = "#666666",
      alpha = 0.75,
      size = label_size,
      max.overlaps = Inf,
      min.segment.length = 0,
      box.padding = 0.5
    )
  } else {
    # Fallback to base geom_text if ggrepel is disabled or unavailable
    scatter <- scatter + geom_text(
      aes(label = IID),
      colour = scales::alpha("#666666", 0.75),
      size = label_size,
      hjust = 0,
      vjust = 0,
      nudge_x = 0.01,
      nudge_y = 0.01
    )
  }
}

ragg::agg_png(file.path(output_dir, "pca_PC1_PC2.png"), width = 1600, height = 1200, res = 200)
print(scatter)
dev.off()

scree <- data.frame(
  PC = seq_along(eigenvalues),
  Proportion = variance_prop
)

scree_plot <- ggplot(scree, aes(x = PC, y = Proportion)) +
  geom_line(colour = "#33a02c", linewidth = 1) +
  geom_point(colour = "#33a02c", size = 2) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  theme_minimal(base_size = 14) +
  labs(
    title = "PCA Scree Plot",
    x = "Principal Component",
    y = "Explained Variance"
  )

ragg::agg_png(file.path(output_dir, "pca_scree.png"), width = 1600, height = 1200, res = 200)
print(scree_plot)
dev.off()

message(sprintf("PCA plots written to %s", normalizePath(output_dir)))
