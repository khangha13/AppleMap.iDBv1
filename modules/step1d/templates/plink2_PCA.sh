#!/bin/bash
# =============================================================================
# Step 1D - PCA Helper
# -----------------------------------------------------------------------------
# Converts combined VCF to PLINK2 format, applies light-weight QC, 
# (optionally) removes close relatives, runs PCA, and generates 
# publication-ready plots.
#
# NOTE: This script expects combined_for_pca.vcf.gz to already exist in the
# VCF source directory. Use prepare_combined_for_pca.sh to create it, or
# the master_vcf_analysis.sh --PCA workflow will auto-prepare it.
# =============================================================================
set -Eeuo pipefail

usage() {
    cat <<'EOF'
Usage: plink2_PCA.sh <vcf_dir> <rscripts_dir> [plink2_bin] [bcftools_bin] \
                     [remove_relatives] [show_labels] [label_size] [use_ggrepel] \
                     [merged_vcf_pattern] [duplicate_mode] [duplicate_king_threshold] \
                     [run_mode]

Arguments:
  vcf_dir              Directory containing combined_for_pca.vcf.gz
  rscripts_dir         Directory containing PCA_plot.R

Optional Arguments (with defaults):
  plink2_bin           Path to plink2 binary (default: plink2)
  bcftools_bin         Path to bcftools binary (default: bcftools)
  remove_relatives     Remove close relatives before PCA (default: false)
  show_labels          Show sample labels on plots (default: true)
  label_size           Plot label size (default: 3)
  use_ggrepel          Use ggrepel for non-overlapping labels (default: true)
  merged_vcf_pattern   Pattern for merged VCF detection (default: *merged*.vcf.gz,*merge*.vcf.gz)
  duplicate_mode       Duplicate handling: off|flag|remove (default: flag)
  duplicate_king_threshold  KING threshold for duplicates (default: 0.45)
  run_mode             Run mode: pca|duplicate (default: pca)

Outputs:
  pca.eigenvec              Eigenvectors (sample coordinates)
  pca.eigenval              Eigenvalues (variance explained)
  pca_plot_PC1_PC2.png      PCA scatter plot
  king_duplicate_pairs.tsv  Duplicate pairs (if detected)
  king_duplicate_samples.tsv Flagged samples (if detected)

Requirements:
  - combined_for_pca.vcf.gz must exist in vcf_dir
  - PLINK2 and bcftools available
  - R with ggplot2, data.table, ragg, scales

Examples:
  # Basic PCA
  plink2_PCA.sh /data/vcfs /path/to/Rscripts

  # With relative removal
  plink2_PCA.sh /data/vcfs /path/to/Rscripts plink2 bcftools true

  # Duplicate detection only
  plink2_PCA.sh /data/vcfs /path/to/Rscripts plink2 bcftools false true 3 true "" flag 0.45 duplicate
EOF
}

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

# Handle --help flag
if [ "${1:-}" = "--help" ] || [ "${1:-}" = "-h" ]; then
    usage
    exit 0
fi

if [ "$#" -lt 2 ]; then
    log_error "Usage: plink2_PCA.sh <vcf_dir> <rscripts_dir> [plink2_bin] [bcftools_bin] [remove_relatives] [show_labels] [label_size] [use_ggrepel] [merged_vcf_pattern] [duplicate_mode] [duplicate_king_threshold] [run_mode]"
    log_error "Run 'plink2_PCA.sh --help' for detailed usage information"
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
MERGED_PATTERN_RAW="${9:-${STEP1D_PCA_MERGED_PATTERN:-}}"
FORCE_CONCAT_RAW="${STEP1D_PCA_FORCE_CONCAT:-false}"
MERGED_EXCLUDE_CHR_RAW="${STEP1D_PCA_MERGED_EXCLUDE_CHR:-true}"
REUSE_COMBINED_RAW="${STEP1D_PCA_REUSE_COMBINED:-true}"
# PLINK QC thresholds (override via environment at submit time)
PCA_GENO_RAW="${STEP1D_PCA_GENO:-0.05}"
PCA_MIND_RAW="${STEP1D_PCA_MIND:-0.10}"
PCA_MAF_RAW="${STEP1D_PCA_MAF:-0.01}"
DUPLICATE_MODE_RAW="${10:-${STEP1D_DUPLICATE_MODE:-flag}}"
DUPLICATE_THRESHOLD_RAW="${11:-${STEP1D_DUPLICATE_KING_THRESHOLD:-0.45}}"
RUN_MODE_RAW="${12:-pca}"

normalize_bool() {
    local value="$1"
    value="$(printf '%s' "${value}" | tr '[:upper:]' '[:lower:]')"
    case "${value}" in
        true|1|yes|y) echo "true" ;;
        *) echo "false" ;;
    esac
}

normalize_duplicate_mode() {
    local value="$1"
    value="$(printf '%s' "${value}" | tr '[:upper:]' '[:lower:]')"
    case "${value}" in
        ""|off|false|0|no|none) echo "off" ;;
        remove|dedup|drop|exclude) echo "remove" ;;
        flag|true|1|yes|on) echo "flag" ;;
        *) echo "${value}" ;;
    esac
}

trim_pattern() {
    printf '%s' "$1" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//'
}

normalize_run_mode() {
    local value="$1"
    value="$(printf '%s' "${value}" | tr '[:upper:]' '[:lower:]')"
    case "${value}" in
        duplicate|dup|duplicates) echo "duplicate" ;;
        pca|"") echo "pca" ;;
        *) echo "${value}" ;;
    esac
}

REMOVE_RELATIVES="$(normalize_bool "${REMOVE_RELATIVES_RAW}")"
SHOW_LABELS="$(normalize_bool "${SHOW_LABELS_RAW}")"
USE_GGREPEL="$(normalize_bool "${USE_GGREPEL_RAW}")"
FORCE_CONCAT="$(normalize_bool "${FORCE_CONCAT_RAW}")"
MERGED_EXCLUDE_CHR="$(normalize_bool "${MERGED_EXCLUDE_CHR_RAW}")"
REUSE_COMBINED="$(normalize_bool "${REUSE_COMBINED_RAW}")"
DUPLICATE_MODE="$(normalize_duplicate_mode "${DUPLICATE_MODE_RAW}")"
RUN_MODE="$(normalize_run_mode "${RUN_MODE_RAW}")"

normalize_num() {
    local value="$1"
    local fallback="$2"
    if [[ "${value}" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
        printf '%s' "${value}"
    else
        printf '%s' "${fallback}"
    fi
}

PCA_GENO="$(normalize_num "${PCA_GENO_RAW}" "0.05")"
PCA_MIND="$(normalize_num "${PCA_MIND_RAW}" "0.10")"
PCA_MAF="$(normalize_num "${PCA_MAF_RAW}" "0.01")"

if [[ "${DUPLICATE_THRESHOLD_RAW}" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
    DUPLICATE_THRESHOLD="${DUPLICATE_THRESHOLD_RAW}"
else
    DUPLICATE_THRESHOLD="0.45"
fi

case "${DUPLICATE_MODE}" in
    off|flag|remove)
        ;;
    *)
        log_warn "Unknown duplicate mode '${DUPLICATE_MODE_RAW}'; defaulting to 'flag'."
        DUPLICATE_MODE="flag"
        ;;
esac

case "${RUN_MODE}" in
    pca|duplicate)
        ;;
    *)
        log_warn "Unknown run mode '${RUN_MODE_RAW}'; defaulting to 'pca'."
        RUN_MODE="pca"
        ;;
esac

if [ "${RUN_MODE}" = "duplicate" ] && [ "${DUPLICATE_MODE}" = "off" ]; then
    log_info "Duplicate-check mode requested; enabling duplicate detection (flag mode)."
    DUPLICATE_MODE="flag"
fi

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

# Expect combined_for_pca.vcf.gz to already exist
# (prepared by master_vcf_analysis.sh or prepare_combined_for_pca.sh)
COMBINED_VCF="combined_for_pca.vcf.gz"

if [ ! -f "${VCF_SOURCE_DIR}/${COMBINED_VCF}" ]; then
    log_error "Required combined VCF not found: ${VCF_SOURCE_DIR}/${COMBINED_VCF}"
    log_error ""
    log_error "This script expects combined_for_pca.vcf.gz to already exist."
    log_error "Please prepare it first using one of these methods:"
    log_error ""
    log_error "  Option 1 (Recommended): Use the PCA workflow which auto-prepares it:"
    log_error "    bash modules/step1d/bin/run_step1d.sh ${VCF_SOURCE_DIR} --PCA"
    log_error ""
    log_error "  Option 2: Manually prepare the combined VCF:"
    log_error "    bash modules/step1d/bin/prepare_combined_for_pca.sh ${VCF_SOURCE_DIR}"
    log_error "    # Then run this script again"
    log_error ""
    exit 1
fi

# Work in PCA output directory, link combined VCF if needed
if [ "${VCF_SOURCE_DIR}/${COMBINED_VCF}" != "${PWD}/${COMBINED_VCF}" ]; then
    log_info "Linking combined VCF from source directory"
    ln -sf "${VCF_SOURCE_DIR}/${COMBINED_VCF}" "${COMBINED_VCF}"
    
    # Link index files (try both CSI and TBI)
    if [ -f "${VCF_SOURCE_DIR}/${COMBINED_VCF}.csi" ]; then
        ln -sf "${VCF_SOURCE_DIR}/${COMBINED_VCF}.csi" "${COMBINED_VCF}.csi"
    fi
    if [ -f "${VCF_SOURCE_DIR}/${COMBINED_VCF}.tbi" ]; then
        ln -sf "${VCF_SOURCE_DIR}/${COMBINED_VCF}.tbi" "${COMBINED_VCF}.tbi"
    fi
fi

log_info "Using combined VCF for PCA: ${COMBINED_VCF}"

IMPORT_PREFIX="all_chromosomes"
QC_PREFIX="qc"

import_stamp="${IMPORT_PREFIX}.import.params.txt"
qc_stamp="${QC_PREFIX}.qc.params.txt"

log_info "Converting VCF to PLINK2 pgen format"
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
        log_info "Reusing existing ${IMPORT_PREFIX}.{pgen,pvar,psam} (import params unchanged; newer than ${COMBINED_VCF})."
    fi
fi

if [ "${needs_import}" = "true" ]; then
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

log_info "Basic QC (geno=${PCA_GENO}, mind=${PCA_MIND}, maf=${PCA_MAF})"
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
    # Reuse only when params match and qc output is newer than import output
    if diff -q <(printf '%s' "${qc_params}") "${qc_stamp}" >/dev/null 2>&1 && [ "${QC_PREFIX}.pgen" -nt "${IMPORT_PREFIX}.pgen" ]; then
        needs_qc="false"
        log_info "Reusing existing ${QC_PREFIX}.{pgen,pvar,psam} (QC params unchanged; newer than ${IMPORT_PREFIX})."
    fi
fi

if [ "${needs_qc}" = "true" ]; then
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

DUPLICATE_SAMPLES_FILE=""
if [ "${DUPLICATE_MODE}" != "off" ]; then
    log_info "Estimating KING kinship for duplicate detection (threshold ${DUPLICATE_THRESHOLD})"
    KING_PREFIX="king_duplicates"
    "${PLINK2_BIN}" \
        --pfile qc \
        --make-king-table \
        --out "${KING_PREFIX}"

    KING_FILE=""
    for candidate in "${KING_PREFIX}.king" "${KING_PREFIX}.kin0"; do
        if [ -f "${candidate}" ]; then
            KING_FILE="${candidate}"
            break
        fi
    done
    if [ -z "${KING_FILE}" ]; then
        log_error "KING output not found for prefix: ${KING_PREFIX}"
        exit 1
    fi

    DUPLICATE_PAIRS_FILE="king_duplicate_pairs.tsv"
    DUPLICATE_SAMPLES_FILE="king_duplicate_samples.tsv"
    DUPLICATE_SAMPLES_TMP="${DUPLICATE_SAMPLES_FILE}.tmp"

    awk -v thr="${DUPLICATE_THRESHOLD}" \
        -v pairs="${DUPLICATE_PAIRS_FILE}" \
        -v samples="${DUPLICATE_SAMPLES_TMP}" '
        BEGIN {
            FS = "[[:space:]]+";
            OFS = "\t";
        }
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                col = tolower($i);
                sub(/^#/, "", col);
                if (col == "fid1") fid1 = i;
                else if (col == "iid1") iid1 = i;
                else if (col == "fid2") fid2 = i;
                else if (col == "iid2") iid2 = i;
                else if (col == "kinship" || col == "king" || col ~ /kinship/) kin = i;
            }
            if (!fid1) fid1 = 1;
            if (!iid1) iid1 = 2;
            if (!fid2) fid2 = 3;
            if (!iid2) iid2 = 4;
            if (!kin) kin = NF;
            print "FID1", "IID1", "FID2", "IID2", "KINSHIP" > pairs;
            next;
        }
        {
            if ($kin == "" || $kin == "NA") {
                next;
            }
            if ($kin + 0 >= thr) {
                f1 = $fid1;
                i1 = $iid1;
                f2 = $fid2;
                i2 = $iid2;
                printf "%s\t%s\t%s\t%s\t%s\n", f1, i1, f2, i2, $kin >> pairs;
                printf "%s\t%s\n", f1, i1 >> samples;
                printf "%s\t%s\n", f2, i2 >> samples;
            }
        }
    ' "${KING_FILE}"

    if [ -s "${DUPLICATE_SAMPLES_TMP}" ]; then
        LC_ALL=C sort -u "${DUPLICATE_SAMPLES_TMP}" > "${DUPLICATE_SAMPLES_FILE}"
    else
        rm -f "${DUPLICATE_SAMPLES_FILE}" 2>/dev/null || true
    fi
    rm -f "${DUPLICATE_SAMPLES_TMP}" 2>/dev/null || true

    DUPLICATE_PAIR_COUNT=0
    if [ -s "${DUPLICATE_PAIRS_FILE}" ]; then
        DUPLICATE_PAIR_COUNT=$(($(wc -l < "${DUPLICATE_PAIRS_FILE}") - 1))
    fi
    DUPLICATE_SAMPLE_COUNT=0
    if [ -s "${DUPLICATE_SAMPLES_FILE}" ]; then
        DUPLICATE_SAMPLE_COUNT=$(wc -l < "${DUPLICATE_SAMPLES_FILE}")
    fi

    if [ "${DUPLICATE_PAIR_COUNT}" -eq 0 ]; then
        log_info "No duplicate-level kinship pairs detected (threshold ${DUPLICATE_THRESHOLD})"
        rm -f "${DUPLICATE_PAIRS_FILE}" "${DUPLICATE_SAMPLES_FILE}" 2>/dev/null || true
        DUPLICATE_SAMPLES_FILE=""
    else
        log_info "Duplicate detection: ${DUPLICATE_PAIR_COUNT} pair(s), ${DUPLICATE_SAMPLE_COUNT} sample(s) flagged"
    fi

    rm -f "${KING_FILE}" "${KING_PREFIX}.log" 2>/dev/null || true

    if [ "${RUN_MODE}" = "duplicate" ]; then
        if [ -n "${DUPLICATE_PAIRS_FILE:-}" ] && [ -f "${DUPLICATE_PAIRS_FILE}" ]; then
            log_info "Duplicate-check complete. Pairs: ${DUPLICATE_PAIR_COUNT}; samples: ${DUPLICATE_SAMPLE_COUNT}"
        else
            log_info "Duplicate-check complete. No pairs above threshold ${DUPLICATE_THRESHOLD}."
        fi
        exit 0
    fi
fi

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
if [ "${DUPLICATE_MODE}" = "remove" ] && [ -n "${DUPLICATE_SAMPLES_FILE}" ] && [ -s "${DUPLICATE_SAMPLES_FILE}" ]; then
    log_info "Removing duplicate samples for PCA"
    "${PLINK2_BIN}" \
        --pfile "${PCA_INPUT}" \
        --remove "${DUPLICATE_SAMPLES_FILE}" \
        --make-pgen \
        --out qc_pruned_dedup
    PCA_INPUT="qc_pruned_dedup"
fi

if [ "${REMOVE_RELATIVES}" = "true" ]; then
    log_info "Removing close relatives (KING cutoff 0.125)"
    "${PLINK2_BIN}" \
        --pfile "${PCA_INPUT}" \
        --king-cutoff 0.125 \
        --out kin

    if [ ! -s "kin.king.cutoff.in.id" ]; then
        log_warn "No individuals passed the KING cutoff; continuing with all samples."
    else
        NOREL_OUT="${PCA_INPUT}_norel"
        "${PLINK2_BIN}" \
            --pfile "${PCA_INPUT}" \
            --keep kin.king.cutoff.in.id \
            --make-pgen \
            --out "${NOREL_OUT}"
        PCA_INPUT="${NOREL_OUT}"
    fi
fi

if [ ! -f "${PCA_INPUT}.pgen" ]; then
    log_error "PCA input dataset not found (${PCA_INPUT}.pgen)"
    exit 1
fi

log_info "Running PCA (10 components with variant weights)"
"${PLINK2_BIN}" \
    --pfile "${PCA_INPUT}" \
    --pca 10 biallelic-var-wts \
    --out pca

log_info "Rendering PCA plots"
PCA_DUPLICATE_FILE=""
if [ -n "${DUPLICATE_SAMPLES_FILE}" ] && [ -s "${DUPLICATE_SAMPLES_FILE}" ]; then
    PCA_DUPLICATE_FILE="${DUPLICATE_SAMPLES_FILE}"
fi
Rscript "${RSCRIPTS_DIR}/PCA_plot.R" "${PWD}/pca.eigenvec" "${PWD}/pca.eigenval" "${PWD}" "${SHOW_LABELS}" "${LABEL_SIZE}" "${USE_GGREPEL}" "${PCA_DUPLICATE_FILE}"

log_info "Cleaning up large intermediate files"
rm -f all_chromosomes.pgen all_chromosomes.pvar all_chromosomes.psam all_chromosomes.log 2>/dev/null || true
rm -f qc.log qc.nosex qc.pgen qc.pvar qc.psam pruned.log pruned.prune.in pruned.prune.out 2>/dev/null || true
rm -f kin.log kin.king.cutoff.in.id kin.king.cutoff.out.id 2>/dev/null || true

log_info "PCA workflow completed."
