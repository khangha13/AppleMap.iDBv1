#!/usr/bin/env Rscript
#
# Boxplot of SNPs per sample by chromosome (chr_label)
# Input: snps_per_sample_counts.csv(.gz) from snp_counts_per_sample_from_chr_vcfs.sh
#
# RStudio usage:
# - Open this file and click Source, or run line-by-line.
# - Set infile/outfile below (or call from terminal via Rscript).
#
# CLI usage:
#   Rscript plot_snps_per_sample_boxplot_by_chr.R snps_per_sample_counts.csv.gz out.png
#

args <- commandArgs(trailingOnly = TRUE)
infile <- if (length(args) >= 1) args[[1]] else "snps_per_sample_counts.csv.gz"
outfile <- if (length(args) >= 2) args[[2]] else "snps_per_sample_boxplot_by_chr.png"

if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")

library(data.table)
library(ggplot2)

dt <- fread(infile)

required <- c("chr_label", "n_called")
missing_cols <- setdiff(required, names(dt))
if (length(missing_cols) > 0) stop("Missing required columns: ", paste(missing_cols, collapse = ", "))

# Keep only per-chromosome rows
dt <- dt[chr_label != "TOTAL"]
dt[, n_called := as.numeric(n_called)]
dt <- dt[!is.na(n_called)]

# Ensure discrete grouping so the box actually renders
dt[, chr_label := factor(chr_label, levels = unique(chr_label))]

p <- ggplot(dt, aes(x = chr_label, y = n_called, group = chr_label)) +
  geom_boxplot(
    fill = "grey85",
    color = "black",
    width = 0.7,
    outlier.shape = 16,
    outlier.size = 1.2
  ) +
  labs(
    x = "Chromosome (from filename)",
    y = "Called SNP genotypes per sample",
    title = "SNPs per sample by chromosome"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# RStudio: display the plot
print(p)

# Optional: save to file (uncomment)
# ggsave(outfile, p, width = 10, height = 4.5, dpi = 200)

