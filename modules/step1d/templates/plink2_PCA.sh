#!/bin/bash
# =============================================================================
# Step 1D - PCA Helper (v2 -- compute everything, decide later)
# -----------------------------------------------------------------------------
# Converts combined VCF to PLINK2 format, applies light-weight QC, computes
# all-pairwise KING kinship, LD-prunes, and runs PCA on ALL samples.
#
# No duplicate filtering, no relative removal, no plotting.
# The raw KING .kin0 file is preserved for downstream Parquet export.
#
# Positional arguments (5):
#   1. vcf_dir        Directory containing combined_for_pca.vcf.gz
#   2. rscripts_dir   Directory containing R scripts (unused here, kept for API)
#   3. plink2_bin     Path to plink2 binary (default: plink2)
#   4. bcftools_bin   Path to bcftools binary (default: bcftools)
#   5. cache_pca_dir  Output directory for all PCA artefacts
# =============================================================================
set -Eeuo pipefail

usage() {
    cat <<'EOF'
Usage: plink2_PCA.sh <vcf_dir> <rscripts_dir> [plink2_bin] [bcftools_bin] [cache_pca_dir]

Arguments:
  vcf_dir         Directory containing combined_for_pca.vcf.gz
  rscripts_dir    Directory containing R scripts (kept for API compatibility)
  plink2_bin      Path to plink2 binary (default: plink2)
  bcftools_bin    Path to bcftools binary (default: bcftools)
  cache_pca_dir   Output directory for PCA artefacts (default: current directory)

Pipeline steps (each with stamp-based skip check):
  1. Import VCF to PLINK (--snps-only just-acgt --max-alleles 2)
  2. QC filter (--geno 0.05 --mind 0.10 --maf 0.01)
  3. KING table (--make-king-table, preserves raw .kin0)
  4. LD prune (--indep-pairwise 200 50 0.2)
  5. PCA on all pruned samples (--pca 10 biallelic-var-wts)

Outputs:
  all_chromosomes.{pgen,pvar,psam}   Imported dataset
  qc.{pgen,pvar,psam}                QC-filtered dataset
  king_pairwise.kin0                  All pairwise KING kinship values
  qc_pruned.{pgen,pvar,psam}         LD-pruned dataset
  pca.eigenvec                        Eigenvectors (sample coordinates)
  pca.eigenval                        Eigenvalues (variance explained)
EOF
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [ -n "${PIPELINE_ROOT:-}" ]; then
    if [ -d "${PIPELINE_ROOT}" ]; then
        PIPELINE_ROOT="$(cd "${PIPELINE_ROOT}" && pwd)"
    else
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
if ! command -v log_success >/dev/null 2>&1; then
    log_success() { log_info "$1"; }
fi

trap 'log_error "Unexpected failure in plink2_PCA.sh"; exit 1' ERR

if [ "${1:-}" = "--help" ] || [ "${1:-}" = "-h" ]; then
    usage
    exit 0
fi

if [ "$#" -lt 2 ]; then
    log_error "Usage: plink2_PCA.sh <vcf_dir> <rscripts_dir> [plink2_bin] [bcftools_bin] [cache_pca_dir]"
    log_error "Run 'plink2_PCA.sh --help' for detailed usage information"
    exit 1
fi

VCF_SOURCE_DIR="$1"
RSCRIPTS_DIR="$2"
PLINK2_BIN="${3:-plink2}"
BCFTOOLS_BIN="${4:-bcftools}"
CACHE_PCA_DIR="${5:-${PWD}}"

# QC thresholds (overridable via environment)
PCA_GENO="${STEP1D_PCA_GENO:-0.05}"
PCA_MIND="${STEP1D_PCA_MIND:-0.10}"
PCA_MAF="${STEP1D_PCA_MAF:-0.01}"

if [ ! -d "${VCF_SOURCE_DIR}" ]; then
    log_error "VCF directory not found: ${VCF_SOURCE_DIR}"
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

COMBINED_VCF="combined_for_pca.vcf.gz"

if [ ! -f "${VCF_SOURCE_DIR}/${COMBINED_VCF}" ]; then
    log_error "Required combined VCF not found: ${VCF_SOURCE_DIR}/${COMBINED_VCF}"
    log_error "Please prepare it first using prepare_combined_for_pca.sh"
    exit 1
fi

mkdir -p "${CACHE_PCA_DIR}"
cd "${CACHE_PCA_DIR}"

# Link combined VCF if not already in our working directory
if [ "${VCF_SOURCE_DIR}/${COMBINED_VCF}" != "${PWD}/${COMBINED_VCF}" ]; then
    log_info "Linking combined VCF from source directory"
    ln -sf "${VCF_SOURCE_DIR}/${COMBINED_VCF}" "${COMBINED_VCF}"
    if [ -f "${VCF_SOURCE_DIR}/${COMBINED_VCF}.csi" ]; then
        ln -sf "${VCF_SOURCE_DIR}/${COMBINED_VCF}.csi" "${COMBINED_VCF}.csi"
    fi
    if [ -f "${VCF_SOURCE_DIR}/${COMBINED_VCF}.tbi" ]; then
        ln -sf "${VCF_SOURCE_DIR}/${COMBINED_VCF}.tbi" "${COMBINED_VCF}.tbi"
    fi
fi

log_info "Using combined VCF for PCA: ${COMBINED_VCF}"

# =============================================================================
# STEP 1: Import VCF to PLINK2
# =============================================================================
IMPORT_PREFIX="all_chromosomes"
import_stamp="${IMPORT_PREFIX}.import.params.txt"
import_params="$(cat <<EOF
vcf=${COMBINED_VCF}
snps_only=just-acgt
max_alleles=2
set_all_var_ids=@:#:\$r:\$a
new_id_max_allele_len=1000
EOF
)"

needs_import="true"
if [ -f "${IMPORT_PREFIX}.pgen" ] && [ -f "${IMPORT_PREFIX}.pvar" ] && [ -f "${IMPORT_PREFIX}.psam" ] && [ -f "${import_stamp}" ]; then
    if diff -q <(printf '%s' "${import_params}") "${import_stamp}" >/dev/null 2>&1 && [ "${IMPORT_PREFIX}.pgen" -nt "${COMBINED_VCF}" ]; then
        needs_import="false"
        log_info "Reusing existing ${IMPORT_PREFIX}.{pgen,pvar,psam} (params unchanged; newer than ${COMBINED_VCF})."
    fi
fi

if [ "${needs_import}" = "true" ]; then
    log_info "Converting VCF to PLINK2 pgen format"
    "${PLINK2_BIN}" \
        --vcf "${COMBINED_VCF}" \
        --double-id \
        --allow-extra-chr \
        --set-all-var-ids '@:#:$r:$a' \
        --snps-only just-acgt \
        --max-alleles 2 \
        --new-id-max-allele-len 1000 \
        --make-pgen \
        --out "${IMPORT_PREFIX}"
    printf '%s' "${import_params}" > "${import_stamp}"
fi

# =============================================================================
# STEP 2: QC filtering
# =============================================================================
QC_PREFIX="qc"
qc_stamp="${QC_PREFIX}.qc.params.txt"
qc_params="$(cat <<EOF
pfile=${IMPORT_PREFIX}
geno=${PCA_GENO}
mind=${PCA_MIND}
maf=${PCA_MAF}
max_alleles=2
EOF
)"

needs_qc="true"
if [ -f "${QC_PREFIX}.pgen" ] && [ -f "${QC_PREFIX}.pvar" ] && [ -f "${QC_PREFIX}.psam" ] && [ -f "${qc_stamp}" ]; then
    if diff -q <(printf '%s' "${qc_params}") "${qc_stamp}" >/dev/null 2>&1 && [ "${QC_PREFIX}.pgen" -nt "${IMPORT_PREFIX}.pgen" ]; then
        needs_qc="false"
        log_info "Reusing existing ${QC_PREFIX}.{pgen,pvar,psam} (QC params unchanged)."
    fi
fi

if [ "${needs_qc}" = "true" ]; then
    log_info "QC filtering (geno=${PCA_GENO}, mind=${PCA_MIND}, maf=${PCA_MAF})"
    "${PLINK2_BIN}" \
        --pfile "${IMPORT_PREFIX}" \
        --geno "${PCA_GENO}" \
        --mind "${PCA_MIND}" \
        --maf "${PCA_MAF}" \
        --max-alleles 2 \
        --make-pgen \
        --out "${QC_PREFIX}"
    printf '%s' "${qc_params}" > "${qc_stamp}"
fi

# =============================================================================
# STEP 3: KING all-pairwise kinship (preserve raw .kin0)
# =============================================================================
KING_PREFIX="king_pairwise"
king_stamp="${KING_PREFIX}.king.params.txt"
king_params="$(cat <<EOF
pfile=${QC_PREFIX}
EOF
)"

needs_king="true"
KING_FILE=""
for candidate in "${KING_PREFIX}.kin0" "${KING_PREFIX}.king"; do
    if [ -f "${candidate}" ]; then
        KING_FILE="${candidate}"
        break
    fi
done

if [ -n "${KING_FILE}" ] && [ -f "${king_stamp}" ]; then
    if diff -q <(printf '%s' "${king_params}") "${king_stamp}" >/dev/null 2>&1 \
       && [ "${KING_FILE}" -nt "${QC_PREFIX}.pgen" ]; then
        needs_king="false"
        log_info "Reusing existing ${KING_FILE} (KING params unchanged)."
    fi
fi

if [ "${needs_king}" = "true" ]; then
    log_info "Computing all-pairwise KING kinship"
    "${PLINK2_BIN}" \
        --pfile "${QC_PREFIX}" \
        --make-king-table \
        --out "${KING_PREFIX}"

    KING_FILE=""
    for candidate in "${KING_PREFIX}.kin0" "${KING_PREFIX}.king"; do
        if [ -f "${candidate}" ]; then
            KING_FILE="${candidate}"
            break
        fi
    done
    if [ -z "${KING_FILE}" ]; then
        log_error "KING output not found for prefix: ${KING_PREFIX}"
        exit 1
    fi

    KING_PAIRS=$(( $(wc -l < "${KING_FILE}") - 1 ))
    log_info "KING computed: ${KING_PAIRS} sample pair(s)"
    printf '%s' "${king_params}" > "${king_stamp}"
fi

# =============================================================================
# STEP 4: LD pruning
# =============================================================================
LD_PRUNED_PREFIX="qc_pruned"
ld_stamp="${LD_PRUNED_PREFIX}.ld.params.txt"
ld_params="$(cat <<EOF
pfile=${QC_PREFIX}
window=200
step=50
r2=0.2
EOF
)"

needs_ld="true"
if [ -f "${LD_PRUNED_PREFIX}.pgen" ] && [ -f "${LD_PRUNED_PREFIX}.pvar" ] && [ -f "${LD_PRUNED_PREFIX}.psam" ] && [ -f "${ld_stamp}" ]; then
    if diff -q <(printf '%s' "${ld_params}") "${ld_stamp}" >/dev/null 2>&1 \
       && [ "${LD_PRUNED_PREFIX}.pgen" -nt "${QC_PREFIX}.pgen" ]; then
        needs_ld="false"
        log_info "Reusing existing ${LD_PRUNED_PREFIX}.{pgen,pvar,psam} (LD params unchanged)."
    fi
fi

if [ "${needs_ld}" = "true" ]; then
    log_info "LD pruning (200 kb window, step 50, r2 0.2)"
    "${PLINK2_BIN}" \
        --pfile "${QC_PREFIX}" \
        --indep-pairwise 200 50 0.2 \
        --out pruned

    "${PLINK2_BIN}" \
        --pfile "${QC_PREFIX}" \
        --extract pruned.prune.in \
        --make-pgen \
        --out "${LD_PRUNED_PREFIX}"
    printf '%s' "${ld_params}" > "${ld_stamp}"
fi

# =============================================================================
# STEP 5: PCA on all pruned samples
# =============================================================================
PCA_INPUT="${LD_PRUNED_PREFIX}"

if [ ! -f "${PCA_INPUT}.pgen" ]; then
    log_error "PCA input dataset not found (${PCA_INPUT}.pgen)"
    exit 1
fi

pca_stamp="pca.pca.params.txt"
pca_params="$(cat <<EOF
pfile=${PCA_INPUT}
pca_components=10
biallelic_var_wts=true
EOF
)"

needs_pca="true"
if [ -f "pca.eigenvec" ] && [ -f "pca.eigenval" ] && [ -f "${pca_stamp}" ]; then
    if diff -q <(printf '%s' "${pca_params}") "${pca_stamp}" >/dev/null 2>&1 \
       && [ "pca.eigenvec" -nt "${PCA_INPUT}.pgen" ]; then
        needs_pca="false"
        log_info "Reusing existing pca.eigenvec/eigenval (PCA params unchanged)."
    fi
fi

if [ "${needs_pca}" = "true" ]; then
    log_info "Running PCA (10 components with variant weights) on all samples"
    "${PLINK2_BIN}" \
        --pfile "${PCA_INPUT}" \
        --pca 10 biallelic-var-wts \
        --out pca
    printf '%s' "${pca_params}" > "${pca_stamp}"
fi

# =============================================================================
# Cleanup: remove only log/temp files; keep data files for skip checks
# =============================================================================
log_info "Cleaning up temporary files (preserving data for future skip checks)"
rm -f "${IMPORT_PREFIX}".log 2>/dev/null || true
rm -f "${QC_PREFIX}".log "${QC_PREFIX}".nosex 2>/dev/null || true
rm -f "${KING_PREFIX}".log 2>/dev/null || true
rm -f pruned.log pruned.prune.in pruned.prune.out 2>/dev/null || true
rm -f "${LD_PRUNED_PREFIX}".log 2>/dev/null || true
rm -f pca.log 2>/dev/null || true

log_info "PCA workflow completed."
