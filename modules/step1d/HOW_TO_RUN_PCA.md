# How to Run PCA Analysis on a VCF File

This guide shows you how to run Principal Component Analysis (PCA) on your VCF files using the Step 1D pipeline.

## Quick Start

### Option 1: Interactive Mode (Easiest - Recommended)

This is the simplest way to run PCA. The script will prompt you for the directory containing your VCF files.

```bash
# Navigate to the pipeline directory
cd /path/to/GATK_Pipeline_KH_v1

# Run PCA-only analysis interactively
bash wrappers/interactive/step1d_interactive.sh --PCA

# If you want to remove close relatives before PCA (optional)
bash wrappers/interactive/step1d_interactive.sh --PCA --remove-relatives
```

**What it does:**
- Prompts you for the directory containing your VCF files
- Automatically detects VCF files (expects `Chr00.vcf.gz` through `Chr17.vcf.gz` by default)
- Runs only PCA analysis (skips QC plots)
- Outputs PCA results in `pca_analysis/` directory

### Option 2: SLURM Batch Mode (For HPC Clusters)

If you're on an HPC cluster with SLURM, use this method:

```bash
# Navigate to the pipeline directory
cd /path/to/GATK_Pipeline_KH_v1

# Submit PCA-only job
bash wrappers/sbatch/step1d_submit.sh <dataset_name> <vcf_directory> --PCA

# Example:
bash wrappers/sbatch/step1d_submit.sh my_dataset /path/to/vcfs --PCA
```

**Parameters:**
- `<dataset_name>`: A name for your dataset (used in job naming)
- `<vcf_directory>`: Full path to directory containing your VCF files
- `--PCA`: Skip QC and run only PCA

### Option 3: Direct Script Call

If you want more control, call the script directly:

```bash
cd /path/to/GATK_Pipeline_KH_v1

# Set environment variables
export VCF_DIR="/path/to/your/vcf/directory"
export VCF_INCLUDE_FILENAMES="your_file.vcf.gz"
export WORK_DIR="${VCF_DIR}"
export STEP1D_REMOVE_RELATIVES=true   # optional

# Run the analysis
bash modules/step1d/templates/master_vcf_analysis.sh --PCA
```

## VCF File Requirements

### Default Naming Pattern
The pipeline expects VCF files named:
- `Chr00.vcf.gz`
- `Chr01.vcf.gz`
- `Chr02.vcf.gz`
- ... through `Chr17.vcf.gz`

**Total: 18 chromosomes (Chr00-Chr17)**

### Custom File Names

If your VCF files have different names, you can specify them:

**For interactive mode:**
```bash
bash wrappers/interactive/step1d_interactive.sh \
  --dir=/path/to/vcfs \
  --vcf=Chr00,Chr01,Chr02 \
  --PCA
```

**For direct script:**
```bash
export VCF_INCLUDE_FILENAMES="Chr00.vcf.gz Chr01.vcf.gz Chr02.vcf.gz"
export VCF_DIR="/path/to/vcfs"
bash modules/step1d/templates/master_vcf_analysis.sh --PCA
```

### Custom VCF Pattern

If your files follow a different pattern (e.g., `chr1.vcf.gz`, `chromosome_01.vcf.gz`):

```bash
export VCF_PATTERN="chr%02d.vcf.gz"  # For chr00.vcf.gz, chr01.vcf.gz, etc.
export VCF_DIR="/path/to/vcfs"
bash modules/step1d/templates/master_vcf_analysis.sh --PCA
```

## Output Files

After running PCA, you'll find results in:

```
<vcf_directory>/pca_analysis/
├── pca.eigenvec          # Principal component vectors (sample coordinates)
├── pca.eigenval          # Eigenvalues (variance explained by each PC)
├── pca.eigenvec.var      # Variance explained per PC
├── pca_PC1_PC2.png       # PCA plot (PC1 vs PC2)
├── pca_PC1_PC3.png       # PCA plot (PC1 vs PC3)
├── pca_PC2_PC3.png       # PCA plot (PC2 vs PC3)
└── pca_scree.png         # Scree plot (variance explained)
```

**Default output directory:** `pca_analysis/` (configurable via `STEP1D_PCA_DIR`)

## Options and Flags

### `--PCA`
- **What it does:** Skips QC stages (depth/missingness/quality plots) and runs only PCA (merged VCF preferred).
- **When to use:** When you only need PCA results and don't need quality control plots.
- **Time saved:** Faster than full QC (plots are skipped).

### `--remove-relatives`
- **What it does:** Removes close relatives (KING coefficient > 0.125) before running PCA
- **Requires:** `--PCA` mode
- **When to use:** When you want to remove related individuals that might bias PCA results
- **Example:**
  ```bash
  bash wrappers/interactive/step1d_interactive.sh --PCA --remove-relatives
  ```

### `--beagle`
- **What it does:** Treats VCFs as Beagle-imputed (expects INFO tags: AF, DR2, IMP)
- **When to use:** If your VCF files were imputed using Beagle
- **Note:** In Beagle mode, expects 17 chromosomes (Chr01-Chr17) instead of 18

### `--dry-run`
- **What it does:** Shows what would be done without actually creating files
- **When to use:** To test your configuration before running the full analysis

## Requirements

### Software Dependencies
The pipeline requires:
- `bcftools` - For VCF manipulation
- `plink2` - For PCA computation
- `R` with packages: `ggplot2`, `data.table`, `ragg`, `scales`, `ggrepel` (optional)

### Module Loading (HPC)
The script automatically loads:
- `miniforge/25.3.0-3`
- `bcftools`
- `plink/2.00a3.6-gcc-11.3.0`

If your HPC uses different module names, you may need to adjust the script or preload modules.

## Examples

### Example 1: Basic PCA on Standard VCF Files
```bash
cd /path/to/GATK_Pipeline_KH_v1
bash wrappers/interactive/step1d_interactive.sh --PCA
# When prompted, enter: /scratch/user/myproject/vcfs
```

### Example 2: PCA with Relative Removal
```bash
cd /path/to/GATK_Pipeline_KH_v1
bash wrappers/interactive/step1d_interactive.sh \
  --dir=/scratch/user/myproject/vcfs \
  --PCA \
  --remove-relatives
```

### Example 3: PCA on Beagle-Imputed VCFs
```bash
cd /path/to/GATK_Pipeline_KH_v1
bash wrappers/interactive/step1d_interactive.sh \
  --dir=/scratch/user/myproject/beagle_vcfs \
  --beagle \
  --PCA
```

### Example 4: Custom VCF File Names
```bash
cd /path/to/GATK_Pipeline_KH_v1
bash wrappers/interactive/step1d_interactive.sh \
  --dir=/scratch/user/myproject/vcfs \
  --vcf=apple_chr1,apple_chr2,apple_chr3 \
  --PCA
```

### Example 5: SLURM Batch Submission
```bash
cd /path/to/GATK_Pipeline_KH_v1
bash wrappers/sbatch/step1d_submit.sh \
  apple_genome \
  /scratch/user/myproject/vcfs \
  --PCA
```

## Troubleshooting

### Problem: "No VCF files specified for analysis"
**Solution:** Check that:
1. Your VCF directory path is correct
2. VCF files are named correctly (default: `Chr00.vcf.gz` through `Chr17.vcf.gz`)
3. Files have `.vcf.gz` or `.vcf` extension

### Problem: "plink2 binary not found"
**Solution:**
```bash
# Load plink2 module
module load plink/2.00a3.6-gcc-11.3.0

# Or specify path manually
export PLINK2_BIN=/path/to/plink2
```

### Problem: "bcftools binary not found"
**Solution:**
```bash
# Load bcftools module
module load bcftools

# Or specify path manually
export BCFTOOLS_BIN=/path/to/bcftools
```

### Problem: "Missing required R packages"
**Solution:**
```bash
# Activate your R environment
conda activate rplot  # or your R environment name

# Install missing packages
Rscript -e "install.packages(c('ggplot2', 'data.table', 'ragg', 'scales', 'ggrepel'), repos='https://cloud.r-project.org')"
```

### Problem: VCF files have different naming
**Solution:** Use `--vcf` flag or set `VCF_INCLUDE_FILENAMES`:
```bash
# Interactive mode
bash wrappers/interactive/step1d_interactive.sh \
  --dir=/path/to/vcfs \
  --vcf=file1,file2,file3 \
  --PCA

# Or export variable
export VCF_INCLUDE_FILENAMES="file1.vcf.gz file2.vcf.gz file3.vcf.gz"
```

## Configuration Options

You can customize PCA behavior with environment variables:

```bash
# Output directory (default: pca_analysis)
export STEP1D_PCA_DIR="my_pca_results"

# Show sample labels on plots (default: true)
export STEP1D_PCA_SHOW_LABELS=true

# Label size (default: 1.5)
export STEP1D_PCA_LABEL_SIZE=2.0

# Use ggrepel for better label positioning (default: true)
export STEP1D_PCA_USE_GGREPEL=true
```

## Time Estimates

- **Small dataset** (<50 samples, <5M variants): 1-3 hours
- **Medium dataset** (50-100 samples, 5-10M variants): 3-6 hours
- **Large dataset** (>100 samples, >10M variants): 6-12+ hours

Times are approximate and depend on:
- Number of samples
- Number of variants
- Cluster performance
- Whether relatives are being removed

## Next Steps

After PCA completes:
1. Check `pca_analysis/pca_PC1_PC2.png` for population structure
2. Review `pca_analysis/pca_scree.png` to see how much variance each PC explains
3. Use `pca_analysis/pca.eigenvec` for downstream analyses
4. Check `pca_analysis/pca.eigenval` for variance explained per component

## Getting Help

- Check logs in the working directory for detailed error messages
- Review `modules/step1d/README_FIRST.md` for general pipeline information
- See `modules/step1d/VCF_Analysis_Quick_Reference.md` for quick commands

## Summary

**Easiest way to run PCA:**
```bash
cd /path/to/GATK_Pipeline_KH_v1
bash wrappers/interactive/step1d_interactive.sh --PCA
```

Enter your VCF directory when prompted, and the pipeline will handle the rest!
