#!/bin/bash
# =============================================================================
# Example: Direct Script Call with Custom VCF Filename and Directory
# =============================================================================
# This shows how to call master_vcf_analysis.sh directly with:
# - Custom VCF filename (not standard Chr00-Chr17 pattern)
# - Custom VCF directory location
# - PCA-only mode
# =============================================================================

# Set the pipeline root (adjust to your actual path)
PIPELINE_ROOT="/path/to/GATK_Pipeline_KH_v1"
# Or let it auto-detect:
# PIPELINE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

# =============================================================================
# REQUIRED: Set your VCF directory and filename
# =============================================================================

# Your custom VCF directory (full path)
export VCF_DIR="/custom/path/to/your/vcf/directory"

# Your VCF filename(s) - space-separated if multiple files
# Option 1: Single VCF file
export VCF_INCLUDE_FILENAMES="my_custom_file.vcf.gz"

# Option 2: Multiple VCF files (space-separated)
# export VCF_INCLUDE_FILENAMES="file1.vcf.gz file2.vcf.gz file3.vcf.gz"

# Option 3: If your files don't have .vcf.gz extension
# export VCF_INCLUDE_FILENAMES="my_file.vcf"

# =============================================================================
# REQUIRED: Set working directory (usually same as VCF_DIR)
# =============================================================================
export WORK_DIR="${VCF_DIR}"

# =============================================================================
# REQUIRED: Enable PCA-only mode
# =============================================================================
export STEP1D_PCA_ONLY=true

# Optional: Remove relatives before PCA
# export STEP1D_REMOVE_RELATIVES=true

# =============================================================================
# OPTIONAL: Set R scripts directory (usually auto-detected)
# =============================================================================
# export R_SCRIPTS_DIR="${PIPELINE_ROOT}/modules/step1d/Rscripts"

# =============================================================================
# OPTIONAL: Customize PCA output directory
# =============================================================================
# export STEP1D_PCA_DIR="my_pca_results"  # Default: pca_analysis

# =============================================================================
# OPTIONAL: PCA plot customization
# =============================================================================
# export STEP1D_PCA_SHOW_LABELS=true      # Default: true
# export STEP1D_PCA_LABEL_SIZE=1.5        # Default: 1.5
# export STEP1D_PCA_USE_GGREPEL=true       # Default: true

# =============================================================================
# OPTIONAL: For Beagle-imputed VCFs
# =============================================================================
# export BEAGLE_MODE=true

# =============================================================================
# Run the script
# =============================================================================
MASTER_SCRIPT="${PIPELINE_ROOT}/modules/step1d/templates/master_vcf_analysis.sh"

if [ ! -f "${MASTER_SCRIPT}" ]; then
    echo "Error: Script not found at ${MASTER_SCRIPT}" >&2
    echo "Please set PIPELINE_ROOT to the correct path" >&2
    exit 1
fi

# Check if VCF file exists
if [ ! -f "${VCF_DIR}/${VCF_INCLUDE_FILENAMES}" ]; then
    echo "Error: VCF file not found: ${VCF_DIR}/${VCF_INCLUDE_FILENAMES}" >&2
    exit 1
fi

echo "Running PCA analysis..."
echo "VCF Directory: ${VCF_DIR}"
echo "VCF File(s): ${VCF_INCLUDE_FILENAMES}"
echo "Output will be in: ${WORK_DIR}/pca_analysis/"
echo ""

# Run with --pca-only flag
bash "${MASTER_SCRIPT}" --pca-only

echo ""
echo "PCA analysis complete!"
echo "Check results in: ${WORK_DIR}/pca_analysis/"
