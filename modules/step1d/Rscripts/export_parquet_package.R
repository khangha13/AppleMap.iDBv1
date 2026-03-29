#!/usr/bin/env Rscript
# =============================================================================
# export_parquet_package.R
# =============================================================================
# Reads Step1D cache outputs (site metrics TSV, PLINK2 PCA files, KING .kin0)
# and writes a self-contained Parquet report package + manifest.json.
#
# Usage:
#   Rscript export_parquet_package.R \
#     --cache-dir /path/to/cache \
#     --package-dir /path/to/output_package \
#     --dataset-name MyDataset \
#     --beagle false \
#     --compression snappy
#
# Dependencies: data.table, arrow, jsonlite
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
  library(jsonlite)
})

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  opts <- list(
    cache_dir    = NULL,
    package_dir  = NULL,
    dataset_name = "unknown",
    beagle       = FALSE,
    compression  = "snappy"
  )
  i <- 1L
  while (i <= length(args)) {
    key <- args[i]
    val <- if (i < length(args)) args[i + 1L] else NULL
    switch(key,
      "--cache-dir"    = { opts$cache_dir    <- val; i <- i + 2L },
      "--package-dir"  = { opts$package_dir  <- val; i <- i + 2L },
      "--dataset-name" = { opts$dataset_name <- val; i <- i + 2L },
      "--beagle"       = { opts$beagle       <- tolower(val) %in% c("true", "1", "yes"); i <- i + 2L },
      "--compression"  = { opts$compression  <- val; i <- i + 2L },
      { stop(sprintf("Unknown argument: %s", key)) }
    )
  }
  if (is.null(opts$cache_dir))   stop("--cache-dir is required")
  if (is.null(opts$package_dir)) stop("--package-dir is required")
  opts
}

opts <- parse_args(args)
message("[export_parquet] Cache:   ", opts$cache_dir)
message("[export_parquet] Package: ", opts$package_dir)
message("[export_parquet] Dataset: ", opts$dataset_name)
message("[export_parquet] Beagle:  ", opts$beagle)
message("[export_parquet] Compression: ", opts$compression)

pca_dir <- file.path(opts$cache_dir, "pca_analysis")

# Create package directories
dir.create(file.path(opts$package_dir, "qc"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(opts$package_dir, "pca"), recursive = TRUE, showWarnings = FALSE)

counts <- list()
sections_present <- character()

# ---------------------------------------------------------------------------
# 1. QC site metrics -> Hive-partitioned Parquet by chrom
# ---------------------------------------------------------------------------
metrics_path <- file.path(opts$cache_dir, "variant_site_metrics.tsv")

if (file.exists(metrics_path)) {
  message("[export_parquet] Reading site metrics: ", metrics_path)
  dt <- fread(metrics_path, sep = "\t", header = TRUE, na.strings = c("NA", ".", ""))

  # Rename columns to lower-snake-case
  old_names <- names(dt)
  new_names <- gsub("([A-Z])", "_\\L\\1", old_names, perl = TRUE)
  new_names <- gsub("^_", "", new_names)
  new_names <- gsub("__+", "_", new_names)
  setnames(dt, old_names, new_names)

  # Ensure unified schema: add missing columns as NA
  standard_cols <- c("chrom", "pos", "qual", "qd", "ac", "af",
                     "inbreeding_coeff", "excess_het", "mq", "mean_depth",
                     "call_rate", "missing_rate", "heterozygous_rate",
                     "dp_non_missing", "called_genotypes", "missing_genotypes",
                     "total_genotypes", "heterozygous_count")
  beagle_cols <- c("chrom", "pos", "qual", "af", "dr2", "imp",
                   "call_rate", "missing_rate", "heterozygous_rate",
                   "called_genotypes", "missing_genotypes",
                   "total_genotypes", "heterozygous_count")
  all_cols <- unique(c(standard_cols, beagle_cols))

  for (col in all_cols) {
    if (!col %in% names(dt)) {
      dt[, (col) := NA]
    }
  }

  # Ensure chrom is character for partitioning
  dt[, chrom := as.character(chrom)]

  chroms <- sort(unique(dt$chrom))
  counts$site_count <- nrow(dt)
  counts$chromosomes <- chroms

  qc_out <- file.path(opts$package_dir, "qc", "site_metrics")
  message("[export_parquet] Writing site metrics (", nrow(dt), " rows, ",
          length(chroms), " chromosomes)")
  write_dataset(dt, qc_out, format = "parquet",
                partitioning = "chrom",
                existing_data_behavior = "overwrite")
  sections_present <- c(sections_present, "qc_site_metrics")
  message("[export_parquet] Site metrics written to ", qc_out)
} else {
  message("[export_parquet] WARNING: Site metrics not found at ", metrics_path)
}

# ---------------------------------------------------------------------------
# 2. PCA scores -> pca/scores.parquet
# ---------------------------------------------------------------------------
eigenvec_path <- file.path(pca_dir, "pca.eigenvec")

if (file.exists(eigenvec_path)) {
  message("[export_parquet] Reading eigenvectors: ", eigenvec_path)
  ev <- fread(eigenvec_path, header = TRUE)

  # PLINK2 eigenvec has #FID IID PC1 PC2 ... PC10
  old_names <- names(ev)
  old_names <- gsub("^#", "", old_names)
  new_names <- tolower(old_names)
  new_names[new_names == "fid"] <- "fid"
  new_names[new_names == "iid"] <- "sample_id"
  setnames(ev, names(ev), new_names)

  counts$sample_count_pca <- nrow(ev)
  counts$pca_components <- sum(grepl("^pc[0-9]+$", names(ev)))

  scores_out <- file.path(opts$package_dir, "pca", "scores.parquet")
  write_parquet(ev, scores_out, compression = opts$compression)
  sections_present <- c(sections_present, "pca_scores")
  message("[export_parquet] PCA scores written (", nrow(ev), " samples, ",
          counts$pca_components, " components)")
} else {
  message("[export_parquet] WARNING: Eigenvector file not found at ", eigenvec_path)
}

# ---------------------------------------------------------------------------
# 3. PCA variance -> pca/variance.parquet
# ---------------------------------------------------------------------------
eigenval_path <- file.path(pca_dir, "pca.eigenval")

if (file.exists(eigenval_path)) {
  message("[export_parquet] Reading eigenvalues: ", eigenval_path)
  eigenvals <- scan(eigenval_path, what = numeric(), quiet = TRUE)
  total_var <- sum(eigenvals)
  var_dt <- data.table(
    component     = seq_along(eigenvals),
    eigenvalue    = eigenvals,
    variance_explained = eigenvals / total_var,
    cumulative_variance_explained = cumsum(eigenvals) / total_var
  )

  variance_out <- file.path(opts$package_dir, "pca", "variance.parquet")
  write_parquet(var_dt, variance_out, compression = opts$compression)
  sections_present <- c(sections_present, "pca_variance")
  message("[export_parquet] Variance explained written (", nrow(var_dt), " components)")
} else {
  message("[export_parquet] WARNING: Eigenvalue file not found at ", eigenval_path)
}

# ---------------------------------------------------------------------------
# 4. Sample annotations -> pca/sample_annotations.parquet
# ---------------------------------------------------------------------------
import_psam <- file.path(pca_dir, "all_chromosomes.psam")
qc_psam <- file.path(pca_dir, "qc.psam")

if (file.exists(import_psam)) {
  message("[export_parquet] Building sample annotations from ", import_psam)
  all_samples <- fread(import_psam, header = TRUE)
  old_names <- gsub("^#", "", names(all_samples))
  setnames(all_samples, names(all_samples), tolower(old_names))

  if ("iid" %in% names(all_samples)) {
    setnames(all_samples, "iid", "sample_id")
  }
  all_samples[, passed_qc := TRUE]
  all_samples[, qc_removal_reason := NA_character_]

  counts$sample_count_import <- nrow(all_samples)

  # Mark samples removed by QC
  if (file.exists(qc_psam)) {
    qc_samples <- fread(qc_psam, header = TRUE)
    qc_names <- gsub("^#", "", names(qc_samples))
    setnames(qc_samples, names(qc_samples), tolower(qc_names))
    if ("iid" %in% names(qc_samples)) {
      setnames(qc_samples, "iid", "sample_id")
    }
    qc_ids <- qc_samples$sample_id
    all_samples[!sample_id %in% qc_ids, `:=`(passed_qc = FALSE, qc_removal_reason = "qc_filter")]
    counts$sample_count_qc <- length(qc_ids)
  }

  # Check for mindrem removal list
  mindrem_path <- file.path(pca_dir, "qc.mindrem.id")
  if (file.exists(mindrem_path)) {
    mindrem <- fread(mindrem_path, header = FALSE)
    if (ncol(mindrem) >= 2) {
      mindrem_ids <- mindrem[[2]]
    } else {
      mindrem_ids <- mindrem[[1]]
    }
    all_samples[sample_id %in% mindrem_ids, `:=`(passed_qc = FALSE, qc_removal_reason = "mind_filter")]
  }

  annot_out <- file.path(opts$package_dir, "pca", "sample_annotations.parquet")
  write_parquet(all_samples, annot_out, compression = opts$compression)
  sections_present <- c(sections_present, "pca_sample_annotations")
  message("[export_parquet] Sample annotations written (", nrow(all_samples), " samples, ",
          sum(!all_samples$passed_qc), " removed by QC)")
} else {
  message("[export_parquet] WARNING: Import .psam not found at ", import_psam)
}

# ---------------------------------------------------------------------------
# 5. King pairwise (all pairs) -> pca/king_pairwise.parquet
# ---------------------------------------------------------------------------
king_kin0 <- NULL
for (candidate in c(file.path(pca_dir, "king_pairwise.kin0"),
                    file.path(pca_dir, "king_pairwise.king"))) {
  if (file.exists(candidate)) {
    king_kin0 <- candidate
    break
  }
}

if (!is.null(king_kin0)) {
  message("[export_parquet] Reading KING table: ", king_kin0)
  king_dt <- fread(king_kin0, header = TRUE)

  # Rename columns to lower-snake-case
  old_names <- gsub("^#", "", names(king_dt))
  new_names <- tolower(old_names)
  setnames(king_dt, names(king_dt), new_names)

  counts$king_pair_count <- nrow(king_dt)

  king_out <- file.path(opts$package_dir, "pca", "king_pairwise.parquet")
  write_parquet(king_dt, king_out, compression = opts$compression)
  sections_present <- c(sections_present, "pca_king_pairwise")
  message("[export_parquet] KING pairwise written (", nrow(king_dt), " pairs)")
} else {
  message("[export_parquet] WARNING: KING .kin0 not found in ", pca_dir)
}

# ---------------------------------------------------------------------------
# 6. Gather counts from PLINK intermediate files
# ---------------------------------------------------------------------------
# Variant counts from .pvar files
for (prefix in c("all_chromosomes", "qc", "qc_pruned")) {
  pvar_path <- file.path(pca_dir, paste0(prefix, ".pvar"))
  if (file.exists(pvar_path)) {
    n_variants <- length(readLines(pvar_path)) - 1L
    key <- paste0("variant_count_", prefix)
    key <- gsub("all_chromosomes", "import", key)
    key <- gsub("qc_pruned", "ld_pruned", key)
    counts[[key]] <- n_variants
  }
}

# ---------------------------------------------------------------------------
# 7. Manifest
# ---------------------------------------------------------------------------
message("[export_parquet] Writing manifest.json")

# QC params from stamp files
qc_params <- list(geno = 0.05, mind = 0.10, maf = 0.01)
qc_stamp <- file.path(pca_dir, "qc.qc.params.txt")
if (file.exists(qc_stamp)) {
  stamp_lines <- readLines(qc_stamp)
  for (line in stamp_lines) {
    parts <- strsplit(line, "=")[[1]]
    if (length(parts) == 2) {
      key <- trimws(parts[1])
      val <- as.numeric(trimws(parts[2]))
      if (!is.na(val) && key %in% c("geno", "mind", "maf")) {
        qc_params[[key]] <- val
      }
    }
  }
}

ld_params <- list(window = 200L, step = 50L, r2 = 0.2)
ld_stamp <- file.path(pca_dir, "qc_pruned.ld.params.txt")
if (file.exists(ld_stamp)) {
  stamp_lines <- readLines(ld_stamp)
  for (line in stamp_lines) {
    parts <- strsplit(line, "=")[[1]]
    if (length(parts) == 2) {
      key <- trimws(parts[1])
      val <- trimws(parts[2])
      if (key %in% c("window", "step")) {
        ld_params[[key]] <- as.integer(val)
      } else if (key == "r2") {
        ld_params[[key]] <- as.numeric(val)
      }
    }
  }
}

manifest <- list(
  bundle_version = "2.0",
  dataset_name   = opts$dataset_name,
  generated_at   = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  input_mode     = ifelse(opts$beagle, "beagle", "standard"),
  chromosomes    = if (exists("chroms")) chroms else character(),
  pca_components = if (!is.null(counts$pca_components)) counts$pca_components else 10L,
  qc_params      = qc_params,
  ld_prune_params = ld_params,
  sections_present = sections_present,
  counts         = counts
)

manifest_path <- file.path(opts$package_dir, "manifest.json")
write_json(manifest, manifest_path, pretty = TRUE, auto_unbox = TRUE)
message("[export_parquet] Manifest written: ", manifest_path)

message("[export_parquet] Export complete. Package: ", opts$package_dir)
