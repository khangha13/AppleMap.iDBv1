#!/usr/bin/env Rscript
# Plot Individual Chromosomes from merged_coverage.bed.gz
# Creates one PNG file per chromosome at 600 dPI

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

# =============================================================================
# CONFIGURATION
# =============================================================================

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  cat("Usage: Rscript plot_individual_chromosomes.R <input.bed.gz> [output_prefix]\n")
  cat("Example: Rscript plot_individual_chromosomes.R merged_coverage.bed.gz chromosomes\n")
  quit(status = 1)
}

input_file <- args[1]
output_prefix <- if (length(args) >= 2) args[2] else "chromosome"

cat("============================================================\n")
cat("Plotting Individual Chromosomes\n")
cat("============================================================\n\n")

cat("Reading coverage data from:", input_file, "\n")
coverage_data <- read.table(input_file, header = FALSE, 
                           col.names = c('chr', 'start', 'end', 'coverage'))

cat("Total data points:", nrow(coverage_data), "\n")

# Filter for chromosomes Chr00-Chr17 (18 chromosomes)
chr_list <- paste0('Chr', sprintf('%02d', 0:17))
coverage_data <- coverage_data[coverage_data$chr %in% chr_list, ]

cat("Data points after filtering:", nrow(coverage_data), "\n\n")

# Convert chromosome to factor for proper ordering
coverage_data$chr <- factor(coverage_data$chr, levels = chr_list)

# Convert position to MegaBasepairs
coverage_data$pos_mb <- coverage_data$start / 1000000

# =============================================================================
# CREATE OUTPUT DIRECTORY
# =============================================================================

output_dir <- paste0(output_prefix, "_plots")
dir.create(output_dir, showWarnings = FALSE)
cat("Output directory:", output_dir, "\n\n")

# =============================================================================
# PLOT EACH CHROMOSOME SEPARATELY
# =============================================================================

cat("Creating individual chromosome plots (600 DPI)...\n\n")

plot_count <- 0

for (chr in chr_list) {
  chr_data <- coverage_data[coverage_data$chr == chr, ]
  
  if (nrow(chr_data) > 0) {
    # Create plot with line thickness 0.3
    p <- ggplot(chr_data, aes(x = pos_mb, y = coverage)) +
      geom_line(color = 'steelblue', linewidth = 0.3, alpha = 0.8) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(size = 12, hjust = 0.5),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        plot.margin = unit(c(0.6, 0.6, 0.6, 0.6), 'cm'),
        panel.grid.minor = element_blank()
      ) +
      labs(
        title = chr,
        subtitle = paste(format(nrow(chr_data), big.mark = ","), "data points"),
        x = 'Position (Mbp)',
        y = 'Coverage (X)'
      ) +
      ylim(0, 200)
    
    # Save individual file at 600 DPI
    output_file <- file.path(output_dir, paste0(chr, ".png"))
    ggsave(output_file, p, width = 10, height = 6, dpi = 600, bg = "white")
    
    plot_count <- plot_count + 1
    cat(sprintf("%2d) Saved: %s.png (%s data points)\n", 
                plot_count, chr, format(nrow(chr_data), big.mark = ",")))
  }
}

cat("\n============================================================\n")
cat("✅ PLOTTING COMPLETE\n")
cat("============================================================\n\n")
cat("Output directory:", output_dir, "\n")
cat("Total plots created:", plot_count, "\n")
cat("Resolution: 600 DPI\n")
cat("Format: Individual PNG files (one per chromosome)\n")
cat("============================================================\n\n")