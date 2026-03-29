# Step1D Output Review (v2 -- report-data)

**Updated:** February 2026
**Purpose:** Document SLURM logging and all output files produced by Step1D

---

## 1. SLURM Output Files (.out and .err)

### Location

SLURM job logs are placed in a dataset-specific subdirectory:

```
${LOG_BASE_PATH}/${dataset_name}/slurm/
```

### File Naming Pattern

```
1D_<JOBID>_<ARRAY_ID>_<JOB_NAME>_<UNIQUE>.output
1D_<JOBID>_<ARRAY_ID>_<JOB_NAME>_<UNIQUE>.error
```

### SLURM Script Location

```
${PIPELINE_SLURM_SCRIPT_DIR}/Apple_GATK_1D_${dataset_name}_<timestamp>.sh
```

---

## 2. Output Directory Structure

Step1D v2 produces two main output directories:

### Cache Directory (intermediate, HPC scratch)

```
${DATASET}_step1d_cache/
  Chr01_snps.vcf.gz ... Chr17_snps.vcf.gz      # Filtered biallelic SNP VCFs
  site_metrics_per_chromosome/                   # Per-chrom metric TSVs
    Chr01_metrics.tsv
    ...
    Chr17_metrics.tsv
  variant_site_metrics.tsv                       # Combined site metrics
  SNP_site_meanDP.tsv                            # Mean depth (standard mode only)
  combined_for_pca.vcf.gz                        # Combined VCF for PCA
  combined_for_pca.vcf.gz.csi                    # Index
  combined_for_pca.stats.txt                     # bcftools stats
  pca_analysis/                                  # PLINK2 intermediates
    all_chromosomes.{pgen,pvar,psam}             # Imported dataset
    all_chromosomes.import.params.txt            # Skip-check stamp
    qc.{pgen,pvar,psam}                          # QC-filtered dataset
    qc.qc.params.txt                             # Skip-check stamp
    king_pairwise.kin0                           # ALL pairwise KING values
    king_pairwise.king.params.txt                # Skip-check stamp
    qc_pruned.{pgen,pvar,psam}                   # LD-pruned dataset
    qc_pruned.ld.params.txt                      # Skip-check stamp
    pca.eigenvec                                 # Sample coordinates (PC1-10)
    pca.eigenval                                 # Eigenvalues
    pca.pca.params.txt                           # Skip-check stamp
```

### Report Package (final output, transfer to local)

```
${DATASET}_step1d_report_package/
  manifest.json                                  # Bundle metadata + counts
  qc/
    site_metrics/                                # Hive-partitioned by chrom
      chrom=Chr01/part-000.parquet
      chrom=Chr02/part-000.parquet
      ...
      chrom=Chr17/part-000.parquet
  pca/
    scores.parquet                               # All samples, PC1-PC10
    variance.parquet                             # 10 rows, eigenvalue + % explained
    sample_annotations.parquet                   # All imported samples + QC flag
    king_pairwise.parquet                        # ALL N*(N-1)/2 pairs
```

Also produced: `${DATASET}_step1d_report_package.tar.gz`

---

## 3. File Descriptions

### Site Metrics TSV

**File:** `variant_site_metrics.tsv`

Standard mode columns:
CHROM, POS, QUAL, QD, AC, AF, INBREEDING_COEFF, EXCESS_HET, MQ, MEAN_DEPTH,
CALL_RATE, MISSING_RATE, HETEROZYGOUS_RATE, DP_NON_MISSING, CALLED_GENOTYPES,
MISSING_GENOTYPES, TOTAL_GENOTYPES, HETEROZYGOUS_COUNT

Beagle mode columns:
CHROM, POS, QUAL, AF, DR2, IMP, CALL_RATE, MISSING_RATE, HETEROZYGOUS_RATE,
CALLED_GENOTYPES, MISSING_GENOTYPES, TOTAL_GENOTYPES, HETEROZYGOUS_COUNT

### Manifest JSON

Contains:
- `bundle_version`: "2.0"
- `dataset_name`: Dataset identifier
- `generated_at`: ISO 8601 timestamp
- `input_mode`: "standard" or "beagle"
- `chromosomes`: Array of chromosome names
- `pca_components`: Number of PCA components (10)
- `qc_params`: geno, mind, maf thresholds
- `ld_prune_params`: window, step, r2
- `sections_present`: Which Parquet sections were successfully exported
- `counts`: site_count, sample_count_import, sample_count_qc, variant counts, king_pair_count

### Parquet Files

| File | Contents |
|------|----------|
| `qc/site_metrics/` | Unified schema (lower-snake-case), missing fields null |
| `pca/scores.parquet` | sample_id, fid, pc1-pc10 |
| `pca/variance.parquet` | component, eigenvalue, variance_explained, cumulative |
| `pca/sample_annotations.parquet` | sample_id, fid, passed_qc, qc_removal_reason |
| `pca/king_pairwise.parquet` | All columns from PLINK2 .kin0 (lower-snake-case) |

---

## 4. Skip-Check System

Each PLINK2 step uses stamp-based skip checks (`.params.txt` files).
A step is skipped if:

1. All output files exist
2. The stamp file exists and matches current parameters
3. Output files are newer than input files

This prevents unnecessary recomputation on reruns.

---

## 5. What Changed from v1

### Removed Outputs

- All plot files (depth, missingness, QD, heterozygosity, AF, call rate, PCA)
- `king_duplicate_pairs.tsv` (threshold filtering is now Quarto-side)
- `king_duplicate_samples.tsv` (same)
- Plot configuration variables

### Added Outputs

- Parquet report package (`manifest.json` + `.parquet` files)
- Tarball of report package
- `king_pairwise.kin0` preserved (all pairs, no filtering)

### Moved Outputs

- All intermediate files now go to `STEP1D_CACHE_DIR` instead of `VCF_DIR`
- `VCF_DIR` is now read-only input

---

## 6. Configuration Variables

| Variable | Default | Description |
|----------|---------|-------------|
| `STEP1D_CACHE_DIR` | `${WORK_DIR}/${DATASET}_step1d_cache` | Intermediate outputs |
| `STEP1D_PACKAGE_DIR` | `${WORK_DIR}/${DATASET}_step1d_report_package` | Final package |
| `STEP1D_PARQUET_COMPRESSION` | `snappy` | Parquet compression |
| `STEP1D_PCA_DIR` | `pca_analysis` | PCA subdirectory name |

---

## Related Files

- **lib/slurm.sh**: SLURM log configuration
- **lib/logging.sh**: Log directory resolution
- **config/pipeline_config.sh**: Central configuration
- **modules/step1d/templates/master_vcf_analysis.sh**: Main orchestrator
- **modules/step1d/templates/plink2_PCA.sh**: PCA pipeline
- **modules/step1d/Rscripts/export_parquet_package.R**: Parquet exporter
