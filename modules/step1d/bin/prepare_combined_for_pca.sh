#!/bin/bash
# =============================================================================
# PREPARE COMBINED VCF FOR PCA (Step 1D Utility)
# =============================================================================
# Detects merged VCF or concatenates per-chromosome VCFs into a single
# combined_for_pca.vcf.gz file with:
#   - Chr00 removed
#   - Contigs standardized to ChrNN format (Chr01..Chr17)
#   - Sorted and indexed (CSI)
#
# Usage:
#   prepare_combined_for_pca.sh <vcf_directory> [--force] [--dry-run] [--help]
#
# Arguments:
#   <vcf_directory>    Directory containing VCF files to process
#
# Options:
#   --force            Overwrite existing combined_for_pca.vcf.gz
#   --dry-run          Show what would be done without creating files
#   --help, -h         Show this help message
#
# Output:
#   <vcf_directory>/combined_for_pca.vcf.gz (+ .csi index)
#   <vcf_directory>/combined_for_pca.stats.txt (bcftools stats output)
#
# NOTE: DATA MANAGEMENT
#   Output files are written into the same directory as the input VCFs
#   (ad-hoc layout). Temporary scratch files use ${TMPDIR:-/tmp}.
#   See the header note in master_vcf_analysis.sh for a full description
#   of the recommended future refactor (separate INPUT/OUTPUT/TEMP paths).
#
# Notes:
#   - If a merged VCF is detected (*merged*.vcf.gz/*merge*.vcf.gz, ignoring
#     names with 'Chr'), it will be cleaned and used.
#   - Otherwise, all *.vcf.gz files in the directory (excluding Chr00) are
#     concatenated.
#   - Contigs are standardized to Chr01..Chr17 format.
#   - Requires bcftools on PATH.
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODULE_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
if [ -n "${PIPELINE_ROOT:-}" ]; then
    if [ -d "${PIPELINE_ROOT}" ]; then
        PIPELINE_ROOT="$(cd "${PIPELINE_ROOT}" && pwd)"
    else
        echo "[prepare_combined_for_pca] ⚠️  Provided PIPELINE_ROOT (${PIPELINE_ROOT}) does not exist; falling back to module-relative path." >&2
        PIPELINE_ROOT="$(cd "${MODULE_DIR}/.." && pwd)"
    fi
else
    PIPELINE_ROOT="$(cd "${MODULE_DIR}/.." && pwd)"
fi
export PIPELINE_ROOT

if [ -f "${PIPELINE_ROOT}/lib/logging.sh" ]; then
    source "${PIPELINE_ROOT}/lib/logging.sh"
    init_logging "prepare_combined_for_pca" "pipeline"
else
    log_info() { echo "[INFO] $*" >&2; }
    log_warn() { echo "[WARN] $*" >&2; }
    log_error() { echo "[ERROR] $*" >&2; }
    log_success() { echo "[SUCCESS] $*" >&2; }
fi

usage() {
    cat <<'EOF'
Usage: prepare_combined_for_pca.sh <vcf_directory> [--force] [--dry-run] [--help]

Arguments:
  <vcf_directory>    Directory containing VCF files to process

Options:
  --force            Overwrite existing combined_for_pca.vcf.gz
  --dry-run          Show what would be done without creating files
  --help, -h         Show this help message

Output:
  <vcf_directory>/combined_for_pca.vcf.gz (+ .csi index)
  <vcf_directory>/combined_for_pca.stats.txt (bcftools stats)

Description:
  Prepares a combined VCF suitable for PCA analysis by:
  1. Detecting a merged VCF or concatenating per-chromosome VCFs
  2. Removing Chr00 records
  3. Standardizing contig names to ChrNN format (Chr01..Chr17)
  4. Sorting and indexing the output
  5. Generating bcftools stats (including total SNP count)

  Selection logic (matches Step1D PCA helper):
  - If a merged VCF is detected (*merged*.vcf.gz/*merge*.vcf.gz, ignoring
    filenames containing 'Chr'), that file is cleaned and used.
  - Otherwise, all *.vcf.gz files (excluding Chr00*.vcf.gz) are concatenated.

Requirements:
  - bcftools (will attempt to load module if not in PATH)
  - tabix (for indexing; included with bcftools)

Examples:
  # Basic usage
  prepare_combined_for_pca.sh /path/to/vcf/directory

  # Overwrite existing output
  prepare_combined_for_pca.sh /path/to/vcf/directory --force

  # Dry run to see what would be done
  prepare_combined_for_pca.sh /path/to/vcf/directory --dry-run
EOF
}

VCF_DIR=""
FORCE=false
DRY_RUN=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --force)
            FORCE=true
            shift
            ;;
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        --help|-h)
            usage
            exit 0
            ;;
        -*)
            log_error "Unknown option: $1"
            usage
            exit 1
            ;;
        *)
            if [ -z "${VCF_DIR}" ]; then
                VCF_DIR="$1"
            else
                log_error "Multiple directory arguments provided"
                usage
                exit 1
            fi
            shift
            ;;
    esac
done

if [ -z "${VCF_DIR}" ]; then
    log_error "VCF directory is required"
    usage
    exit 1
fi

if [ ! -d "${VCF_DIR}" ]; then
    log_error "VCF directory not found: ${VCF_DIR}"
    exit 1
fi

VCF_DIR="$(cd "${VCF_DIR}" && pwd)"
OUTPUT_VCF="${VCF_DIR}/combined_for_pca.vcf.gz"
OUTPUT_STATS="${VCF_DIR}/combined_for_pca.stats.txt"

if [ -f "${OUTPUT_VCF}" ] && [ "${FORCE}" != "true" ]; then
    log_error "Output file already exists: ${OUTPUT_VCF}"
    log_error "Use --force to overwrite or remove the file manually"
    exit 1
fi

# Load bcftools if not available
if ! command -v bcftools >/dev/null 2>&1; then
    if command -v module >/dev/null 2>&1; then
        log_info "bcftools not found, attempting to load module..."
        module purge
        if module load bcftools/1.18-gcc-12.3.0 >/dev/null 2>&1; then
            log_success "Loaded bcftools/1.18-gcc-12.3.0 module"
        else
            log_error "Failed to load bcftools module"
            log_error "Please ensure bcftools is available on PATH"
            exit 1
        fi
    else
        log_error "bcftools not found and module system not available"
        exit 1
    fi
fi

if ! command -v bcftools >/dev/null 2>&1; then
    log_error "bcftools not found in PATH after module load attempt"
    exit 1
fi

BCFTOOLS_BIN="bcftools"

log_info "VCF directory: ${VCF_DIR}"
log_info "Output VCF: ${OUTPUT_VCF}"
log_info "Output stats: ${OUTPUT_STATS}"
if [ "${DRY_RUN}" = "true" ]; then
    log_info "DRY RUN MODE - no files will be created"
fi

# Create temp directory
if [ "${DRY_RUN}" != "true" ]; then
    TEMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/prepare_combined_pca.XXXXXX")"
    trap 'rm -rf "${TEMP_DIR}"' EXIT
    log_info "Temporary directory: ${TEMP_DIR}"
else
    TEMP_DIR="/tmp/dry_run_temp"
    log_info "[dry-run] Would create temporary directory"
fi

# Detect merged VCF (Step1D logic)
MERGED_VCF=""
MERGED_PATTERN="${STEP1D_PCA_MERGED_PATTERN:-*merged*.vcf.gz,*merge*.vcf.gz}"
MERGED_EXCLUDE_CHR="${STEP1D_PCA_MERGED_EXCLUDE_CHR:-true}"

log_info "Searching for merged VCF using pattern: ${MERGED_PATTERN}"
if [ "${MERGED_EXCLUDE_CHR}" = "true" ]; then
    log_info "Will exclude filenames containing 'Chr' from merged VCF detection"
fi

IFS=',' read -r -a MERGED_PATTERNS <<< "${MERGED_PATTERN}"
for pattern in "${MERGED_PATTERNS[@]}"; do
    pattern="$(echo "${pattern}" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')"
    [ -z "${pattern}" ] && continue
    
    mapfile -d '' -t merged_candidates < <(LC_ALL=C find "${VCF_DIR}" -maxdepth 1 \( -type f -o -type l \) -iname "${pattern}" -print0 | LC_ALL=C sort -z) || true
    
    if [ ${#merged_candidates[@]} -gt 0 ]; then
        if [ "${MERGED_EXCLUDE_CHR}" = "true" ]; then
            filtered_candidates=()
            for candidate in "${merged_candidates[@]}"; do
                candidate_base="$(basename "${candidate}")"
                if [[ "${candidate_base}" =~ [Cc][Hh][Rr] ]]; then
                    log_info "Excluding from merged detection (contains 'Chr'): ${candidate_base}"
                    continue
                fi
                filtered_candidates+=("${candidate}")
            done
            merged_candidates=("${filtered_candidates[@]}")
        fi
        
        if [ ${#merged_candidates[@]} -gt 0 ]; then
            MERGED_VCF="${merged_candidates[0]}"
            break
        fi
    fi
done

# Determine input VCFs
declare -a INPUT_VCFS=()
if [ -n "${MERGED_VCF}" ]; then
    log_info "Detected merged VCF: ${MERGED_VCF}"
    INPUT_VCFS=("${MERGED_VCF}")
else
    log_info "No merged VCF detected, will concatenate per-chromosome VCFs"
    
    # Find all *.vcf.gz files, excluding Chr00
    mapfile -d '' -t all_vcfs < <(LC_ALL=C find "${VCF_DIR}" -maxdepth 1 \( -type f -o -type l \) -name "*.vcf.gz" -print0 | LC_ALL=C sort -z) || true
    
    for vcf in "${all_vcfs[@]}"; do
        base="$(basename "${vcf}")"
        # Skip Chr00 files
        if [[ "${base}" =~ ^Chr00 ]]; then
            log_info "Excluding Chr00 file: ${base}"
            continue
        fi
        INPUT_VCFS+=("${vcf}")
    done
    
    if [ ${#INPUT_VCFS[@]} -eq 0 ]; then
        log_error "No VCF files found in ${VCF_DIR}"
        exit 1
    fi
    
    log_info "Found ${#INPUT_VCFS[@]} VCF file(s) to concatenate"
fi

# Create contig rename map (covers multiple input formats -> ChrNN)
# Handles: 1..17, chr1..chr17, Chr1..Chr17, CHR1..CHR17 -> Chr01..Chr17
RENAME_MAP="${TEMP_DIR}/chr_rename.map"
if [ "${DRY_RUN}" != "true" ]; then
    : > "${RENAME_MAP}"
    for num in $(seq 1 17); do
        padded="$(printf "Chr%02d" "${num}")"
        # Numeric format
        printf "%d\t%s\n" "${num}" "${padded}" >> "${RENAME_MAP}"
        # chr format
        printf "chr%d\t%s\n" "${num}" "${padded}" >> "${RENAME_MAP}"
        # Chr non-padded
        printf "Chr%d\t%s\n" "${num}" "${padded}" >> "${RENAME_MAP}"
        # CHR format
        printf "CHR%d\t%s\n" "${num}" "${padded}" >> "${RENAME_MAP}"
    done
    log_info "Created contig rename map: ${RENAME_MAP}"
else
    log_info "[dry-run] Would create contig rename map covering: 1-17, chr1-chr17, Chr1-Chr17, CHR1-CHR17 -> Chr01-Chr17"
fi

# Process VCFs
if [ ${#INPUT_VCFS[@]} -eq 1 ]; then
    # Single file (merged VCF case)
    INPUT_VCF="${INPUT_VCFS[0]}"
    log_info "Processing single VCF: $(basename "${INPUT_VCF}")"
    
    if [ "${DRY_RUN}" = "true" ]; then
        log_info "[dry-run] Would:"
        log_info "[dry-run]   1. Create tabix index on source VCF (workaround for undefined contigs)"
        log_info "[dry-run]   2. Filter to exclude Chr00 records: bcftools view -t ^Chr00"
        log_info "[dry-run]   3. Rename contigs: bcftools annotate --rename-chrs"
        log_info "[dry-run]   4. Sort and write to ${OUTPUT_VCF}"
        log_info "[dry-run]   5. Create CSI index"
        log_info "[dry-run]   6. Generate bcftools stats: ${OUTPUT_STATS}"
    else
        # QUILT2 workaround: create tabix index before rename to handle undefined contigs
        log_info "Creating tabix index on source VCF (workaround for undefined contig headers)..."
        if ! tabix -f -p vcf "${INPUT_VCF}" 2>/dev/null; then
            log_warn "Failed to create tabix index on ${INPUT_VCF}, proceeding anyway"
        fi
        
        FILTERED_VCF="${TEMP_DIR}/filtered.vcf.gz"
        RENAMED_VCF="${TEMP_DIR}/renamed.vcf.gz"
        
        # Step 1: Filter out Chr00
        log_info "Filtering to exclude Chr00 records..."
        if ! "${BCFTOOLS_BIN}" view -t ^Chr00 -Oz -o "${FILTERED_VCF}" "${INPUT_VCF}"; then
            # If Chr00 doesn't exist or filter fails, just copy the file
            log_warn "Failed to filter Chr00 (may not exist), using original file"
            cp "${INPUT_VCF}" "${FILTERED_VCF}"
        fi
        
        # Step 2: Rename contigs
        log_info "Standardizing contig names to ChrNN format..."
        if ! "${BCFTOOLS_BIN}" annotate --rename-chrs "${RENAME_MAP}" -Oz -o "${RENAMED_VCF}" "${FILTERED_VCF}"; then
            log_error "Failed to rename contigs"
            exit 1
        fi
        
        # Step 3: Sort
        log_info "Sorting and writing final VCF..."
        if ! "${BCFTOOLS_BIN}" sort -Oz -o "${OUTPUT_VCF}" "${RENAMED_VCF}"; then
            log_error "Failed to sort VCF"
            exit 1
        fi
        
        # Step 4: Index
        log_info "Creating CSI index..."
        if ! "${BCFTOOLS_BIN}" index -f -c "${OUTPUT_VCF}"; then
            log_error "Failed to create CSI index"
            exit 1
        fi
        
        # Step 5: Generate stats
        log_info "Generating bcftools stats..."
        if ! "${BCFTOOLS_BIN}" stats "${OUTPUT_VCF}" > "${OUTPUT_STATS}"; then
            log_error "Failed to generate bcftools stats"
            exit 1
        fi
        log_success "Created stats file: ${OUTPUT_STATS}"
        
        log_success "Created combined VCF: ${OUTPUT_VCF}"
    fi
else
    # Multiple files (concatenation case)
    log_info "Processing ${#INPUT_VCFS[@]} VCF files for concatenation"
    
    if [ "${DRY_RUN}" = "true" ]; then
        log_info "[dry-run] Would process each VCF:"
        for vcf in "${INPUT_VCFS[@]}"; do
            log_info "[dry-run]   - $(basename "${vcf}")"
        done
        log_info "[dry-run] Steps for each file:"
        log_info "[dry-run]   1. Create tabix index (workaround)"
        log_info "[dry-run]   2. Rename contigs to ChrNN"
        log_info "[dry-run]   3. Sort"
        log_info "[dry-run] Then concatenate all, filter Chr00, create final output, and generate stats"
    else
        # Process each VCF: tabix, rename, sort
        declare -a PROCESSED_VCFS=()
        for i in "${!INPUT_VCFS[@]}"; do
            vcf="${INPUT_VCFS[$i]}"
            base="$(basename "${vcf}" .vcf.gz)"
            log_info "Processing ${vcf} (${i}/${#INPUT_VCFS[@]})..."
            
            # Create tabix index (QUILT2 workaround)
            if ! tabix -f -p vcf "${vcf}" 2>/dev/null; then
                log_warn "Failed to create tabix index on ${vcf}, proceeding anyway"
            fi
            
            renamed="${TEMP_DIR}/${base}.renamed.vcf.gz"
            sorted="${TEMP_DIR}/${base}.sorted.vcf.gz"
            
            # Rename contigs
            if ! "${BCFTOOLS_BIN}" annotate --rename-chrs "${RENAME_MAP}" -Oz -o "${renamed}" "${vcf}"; then
                log_error "Failed to rename contigs in ${vcf}"
                exit 1
            fi
            
            # Sort
            if ! "${BCFTOOLS_BIN}" sort -Oz -o "${sorted}" "${renamed}"; then
                log_error "Failed to sort ${renamed}"
                exit 1
            fi
            
            # Index
            if ! "${BCFTOOLS_BIN}" index -f -c "${sorted}"; then
                log_error "Failed to index ${sorted}"
                exit 1
            fi
            
            PROCESSED_VCFS+=("${sorted}")
        done
        
        # Concatenate all processed VCFs
        log_info "Concatenating ${#PROCESSED_VCFS[@]} processed VCFs..."
        concat_temp="${TEMP_DIR}/concatenated.vcf.gz"
        if ! "${BCFTOOLS_BIN}" concat -Oz -o "${concat_temp}" "${PROCESSED_VCFS[@]}"; then
            log_error "Failed to concatenate VCFs"
            exit 1
        fi
        
        # Filter out Chr00 from concatenated result
        log_info "Filtering concatenated VCF to exclude Chr00..."
        filtered_concat="${TEMP_DIR}/filtered_concat.vcf.gz"
        if ! "${BCFTOOLS_BIN}" view -t ^Chr00 -Oz -o "${filtered_concat}" "${concat_temp}"; then
            log_warn "Failed to filter Chr00 (may not exist), using concatenated file as-is"
            cp "${concat_temp}" "${filtered_concat}"
        fi
        
        # Final sort and output
        log_info "Final sort and writing to ${OUTPUT_VCF}..."
        if ! "${BCFTOOLS_BIN}" sort -Oz -o "${OUTPUT_VCF}" "${filtered_concat}"; then
            log_error "Failed to sort final VCF"
            exit 1
        fi
        
        # Index final output
        log_info "Creating CSI index on final output..."
        if ! "${BCFTOOLS_BIN}" index -f -c "${OUTPUT_VCF}"; then
            log_error "Failed to create CSI index on final output"
            exit 1
        fi
        
        # Generate stats
        log_info "Generating bcftools stats..."
        if ! "${BCFTOOLS_BIN}" stats "${OUTPUT_VCF}" > "${OUTPUT_STATS}"; then
            log_error "Failed to generate bcftools stats"
            exit 1
        fi
        log_success "Created stats file: ${OUTPUT_STATS}"
        
        log_success "Created combined VCF: ${OUTPUT_VCF}"
    fi
fi

# Summary
if [ "${DRY_RUN}" != "true" ]; then
    log_info "Verifying output..."
    
    # Check for Chr01-Chr17 contigs in header
    contigs_found=$("${BCFTOOLS_BIN}" view -h "${OUTPUT_VCF}" | grep -c "^##contig=<ID=Chr[0-9]" || true)
    log_info "Contigs in output: ${contigs_found} Chr* contigs found"
    
    # Check for Chr00
    if "${BCFTOOLS_BIN}" view -h "${OUTPUT_VCF}" | grep -q "^##contig=<ID=Chr00"; then
        log_warn "Chr00 found in output header (may need manual cleanup)"
    else
        log_success "Chr00 successfully excluded from output"
    fi
    
    # Extract key stats from bcftools stats output
    if [ -f "${OUTPUT_STATS}" ]; then
        snp_count=$(grep "^SN" "${OUTPUT_STATS}" | grep "number of SNPs:" | awk '{print $NF}' || echo "0")
        total_records=$(grep "^SN" "${OUTPUT_STATS}" | grep "number of records:" | awk '{print $NF}' || echo "0")
        sample_count=$(grep "^SN" "${OUTPUT_STATS}" | grep "number of samples:" | awk '{print $NF}' || echo "0")
        
        log_info "=== VCF Statistics ==="
        log_info "Total records: ${total_records}"
        log_info "Total SNPs: ${snp_count}"
        log_info "Samples: ${sample_count}"
        log_info "Stats file: ${OUTPUT_STATS}"
    else
        # Fallback if stats file missing
        variant_count=$("${BCFTOOLS_BIN}" view -H "${OUTPUT_VCF}" | wc -l || echo "0")
        log_info "Total variants in combined VCF: ${variant_count}"
    fi
    
    # File size
    file_size=$(du -h "${OUTPUT_VCF}" | cut -f1)
    log_info "Output file size: ${file_size}"
    
    echo ""
    log_success "✓ Combined VCF prepared successfully"
    log_info "Output VCF: ${OUTPUT_VCF}"
    log_info "Index: ${OUTPUT_VCF}.csi"
    log_info "Stats: ${OUTPUT_STATS}"
    echo ""
    log_info "Ready for PCA analysis with Step1D:"
    log_info "  bash modules/step1d/templates/master_vcf_analysis.sh --PCA"
fi

exit 0
