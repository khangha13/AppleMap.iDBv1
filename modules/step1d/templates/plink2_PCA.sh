#!/bin/bash
# =============================================================================
# Step 1D - PCA Helper
# -----------------------------------------------------------------------------
# Combines per-chromosome VCFs, converts them to PLINK2 format, applies
# light-weight QC, (optionally) removes close relatives, runs PCA, and
# generates publication-ready plots.
# =============================================================================
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [ -n "${PIPELINE_ROOT:-}" ]; then
    if [ -d "${PIPELINE_ROOT}" ]; then
        PIPELINE_ROOT="$(cd "${PIPELINE_ROOT}" && pwd)"
    else
        echo "[plink2_PCA] ⚠️  Provided PIPELINE_ROOT (${PIPELINE_ROOT}) does not exist; falling back to template-relative path." >&2
        PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
    fi
else
    PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
fi
export PIPELINE_ROOT
if [ -f "${PIPELINE_ROOT}/lib/logging.sh" ]; then
    source "${PIPELINE_ROOT}/lib/logging.sh"
    init_logging "step1d_pca" "pipeline"
fi
trap 'log_error "Unexpected failure in plink2_PCA.sh"; exit 1' ERR

if [ "$#" -lt 2 ]; then
    log_error "Usage: plink2_PCA.sh <vcf_dir> <rscripts_dir> [plink2_bin] [bcftools_bin] [remove_relatives] [show_labels] [label_size] [use_ggrepel]"
    exit 1
fi

VCF_SOURCE_DIR="$1"
RSCRIPTS_DIR="$2"
PLINK2_BIN="${3:-plink2}"
BCFTOOLS_BIN="${4:-bcftools}"
REMOVE_RELATIVES_RAW="${5:-false}"
SHOW_LABELS_RAW="${6:-true}"
LABEL_SIZE_RAW="${7:-3}"
USE_GGREPEL_RAW="${8:-true}"

normalize_bool() {
    local value="$1"
    value="$(printf '%s' "${value}" | tr '[:upper:]' '[:lower:]')"
    case "${value}" in
        true|1|yes|y) echo "true" ;;
        *) echo "false" ;;
    esac
}

REMOVE_RELATIVES="$(normalize_bool "${REMOVE_RELATIVES_RAW}")"
SHOW_LABELS="$(normalize_bool "${SHOW_LABELS_RAW}")"
USE_GGREPEL="$(normalize_bool "${USE_GGREPEL_RAW}")"

# Validate numeric label size; default to 3 if invalid
if [[ "${LABEL_SIZE_RAW}" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
    LABEL_SIZE="${LABEL_SIZE_RAW}"
else
    LABEL_SIZE="3"
fi

if [ ! -d "${VCF_SOURCE_DIR}" ]; then
    log_error "VCF directory not found: ${VCF_SOURCE_DIR}"
    exit 1
fi

if [ ! -d "${RSCRIPTS_DIR}" ] || [ ! -f "${RSCRIPTS_DIR}/PCA_plot.R" ]; then
    log_error "PCA_plot.R not found in ${RSCRIPTS_DIR}"
    exit 1
fi

if ! command -v "${PLINK2_BIN}" >/dev/null 2>&1; then
    log_error "plink2 binary not found: ${PLINK2_BIN}"
    exit 1
fi

if ! command -v "${BCFTOOLS_BIN}" >/dev/null 2>&1; then
    log_error "bcftools binary not found: ${BCFTOOLS_BIN}"
    exit 1
fi

declare -a FILTERED_VCFS=()
mapfile -d '' -t FILTERED_VCFS < <(LC_ALL=C find "${VCF_SOURCE_DIR}" -maxdepth 1 -type f -name '*_snps.vcf.gz' -print0 | LC_ALL=C sort -z) || true

if [ ${#FILTERED_VCFS[@]} -eq 0 ]; then
    mapfile -d '' -t FILTERED_VCFS < <(LC_ALL=C find "${VCF_SOURCE_DIR}" -maxdepth 1 -type f -name '*.vcf.gz' ! -name '*_snps.vcf.gz' -print0 | LC_ALL=C sort -z) || true
fi

if [ ${#FILTERED_VCFS[@]} -eq 0 ]; then
    mapfile -d '' -t FILTERED_VCFS < <(LC_ALL=C find "${VCF_SOURCE_DIR}" -maxdepth 1 -type f -name '*.vcf' -print0 | LC_ALL=C sort -z) || true
fi

if [ ${#FILTERED_VCFS[@]} -eq 0 ]; then
    log_error "No VCF files found in ${VCF_SOURCE_DIR}"
    exit 1
fi

log_info "Combining ${#FILTERED_VCFS[@]} VCF file(s) for PCA:"
printf '   • %s\n' "${FILTERED_VCFS[@]}"

COMBINED_VCF="combined_for_pca.vcf.gz"
"${BCFTOOLS_BIN}" concat -Oz -o "${COMBINED_VCF}" "${FILTERED_VCFS[@]}"
"${BCFTOOLS_BIN}" index -f "${COMBINED_VCF}"

log_info "Converting VCF to PLINK2 pgen format"
"${PLINK2_BIN}" \
    --vcf "${COMBINED_VCF}" \
    --double-id \
    --allow-extra-chr \
    --set-all-var-ids '@:#:$r:$a' \
    --make-pgen \
    --out all_chromosomes

log_info "Basic QC (geno=0.05, mind=0.10, maf=0.01)"
"${PLINK2_BIN}" \
    --pfile all_chromosomes \
    --geno 0.05 \
    --mind 0.10 \
    --maf 0.01 \
    --max-alleles 2 \
    --make-pgen \
    --out qc

log_info "LD pruning (200 kb window, step 50, r² 0.2)"
"${PLINK2_BIN}" \
    --pfile qc \
    --indep-pairwise 200 50 0.2 \
    --out pruned

"${PLINK2_BIN}" \
    --pfile qc \
    --extract pruned.prune.in \
    --make-pgen \
    --out qc_pruned

PCA_INPUT="qc_pruned"
if [ "${REMOVE_RELATIVES}" = "true" ]; then
    log_info "Removing close relatives (KING cutoff 0.125)"
    "${PLINK2_BIN}" \
        --pfile qc_pruned \
        --king-cutoff 0.125 \
        --out kin

    if [ ! -s "kin.king.cutoff.in.id" ]; then
        log_warn "No individuals passed the KING cutoff; continuing with all samples."
    else
        "${PLINK2_BIN}" \
            --pfile qc_pruned \
            --keep kin.king.cutoff.in.id \
            --make-pgen \
            --out qc_pruned_norel
        PCA_INPUT="qc_pruned_norel"
    fi
fi

if [ ! -f "${PCA_INPUT}.pgen" ]; then
    log_error "PCA input dataset not found (${PCA_INPUT}.pgen)"
    exit 1
fi

log_info "Running PCA (20 components)"
"${PLINK2_BIN}" \
    --pfile "${PCA_INPUT}" \
    --pca 20 \
    --out pca

log_info "Rendering PCA plots"
Rscript "${RSCRIPTS_DIR}/PCA_plot.R" "${PWD}/pca.eigenvec" "${PWD}/pca.eigenval" "${PWD}" "${SHOW_LABELS}" "${LABEL_SIZE}" "${USE_GGREPEL}"

log_info "Cleaning up large intermediate files"
rm -f "${COMBINED_VCF}" "${COMBINED_VCF}.tbi" 2>/dev/null || true
rm -f all_chromosomes.pgen all_chromosomes.pvar all_chromosomes.psam all_chromosomes.log 2>/dev/null || true
rm -f qc.log qc.nosex qc.pgen qc.pvar qc.psam pruned.log pruned.prune.in pruned.prune.out 2>/dev/null || true
rm -f kin.log kin.king.cutoff.in.id kin.king.cutoff.out.id 2>/dev/null || true

log_info "PCA workflow completed."
