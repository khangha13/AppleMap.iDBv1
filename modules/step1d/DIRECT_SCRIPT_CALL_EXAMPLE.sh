#!/bin/bash
# =============================================================================
# Example: Direct Script Call with Custom VCF Filename and Directory
# =============================================================================
# Shows how to call master_vcf_analysis.sh directly (v2 report-data mode).
# The script runs a single workflow: metrics + KING + PCA + Parquet export.
# =============================================================================

# Set the pipeline root (adjust to your actual path)
PIPELINE_ROOT="/path/to/GATK_Pipeline_KH_v1"

# =============================================================================
# REQUIRED: Set your VCF directory and filename
# =============================================================================

export VCF_DIR="/custom/path/to/your/vcf/directory"

# Your VCF filename(s) - space-separated if multiple files
export VCF_INCLUDE_FILENAMES="my_custom_file.vcf.gz"

# =============================================================================
# REQUIRED: Set working directory (usually same as VCF_DIR)
# =============================================================================
export WORK_DIR="${VCF_DIR}"

# =============================================================================
# OPTIONAL: Cache and package directories
# =============================================================================
# export STEP1D_CACHE_DIR="${WORK_DIR}/my_cache"
# export STEP1D_PACKAGE_DIR="${WORK_DIR}/my_report_package"

# =============================================================================
# OPTIONAL: PCA output subdirectory name
# =============================================================================
# export STEP1D_PCA_DIR="my_pca_results"  # Default: pca_analysis

# =============================================================================
# OPTIONAL: For Beagle-imputed VCFs, add --beagle flag below
# =============================================================================

# =============================================================================
# Run the script
# =============================================================================
MASTER_SCRIPT="${PIPELINE_ROOT}/modules/step1d/templates/master_vcf_analysis.sh"

if [ ! -f "${MASTER_SCRIPT}" ]; then
    echo "Error: Script not found at ${MASTER_SCRIPT}" >&2
    echo "Please set PIPELINE_ROOT to the correct path" >&2
    exit 1
fi

echo "Running Step1D report-data pipeline..."
echo "VCF Directory: ${VCF_DIR}"
echo "VCF File(s): ${VCF_INCLUDE_FILENAMES}"
echo ""

bash "${MASTER_SCRIPT}"
# For Beagle mode: bash "${MASTER_SCRIPT}" --beagle
# For dry-run:     bash "${MASTER_SCRIPT}" --dry-run

echo ""
echo "Step1D complete!"
