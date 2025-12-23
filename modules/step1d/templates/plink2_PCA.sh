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
    log_error "Usage: plink2_PCA.sh <vcf_dir> <rscripts_dir> [plink2_bin] [bcftools_bin] [remove_relatives] [show_labels] [label_size] [use_ggrepel] [merged_vcf_pattern] [duplicate_mode] [duplicate_king_threshold] [run_mode]"
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
DUPLICATE_MODE="$(normalize_duplicate_mode "${DUPLICATE_MODE_RAW}")"
RUN_MODE="$(normalize_run_mode "${RUN_MODE_RAW}")"

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

COMBINED_VCF=""
MERGED_VCF=""
if [ -n "${MERGED_PATTERN_RAW}" ]; then
    IFS=',' read -r -a MERGED_PATTERNS <<< "${MERGED_PATTERN_RAW}"
    for pattern in "${MERGED_PATTERNS[@]}"; do
        pattern="$(trim_pattern "${pattern}")"
        [ -z "${pattern}" ] && continue
        mapfile -d '' -t merged_candidates < <(LC_ALL=C find "${VCF_SOURCE_DIR}" -maxdepth 1 -type f -iname "${pattern}" -print0 | LC_ALL=C sort -z) || true
        if [ ${#merged_candidates[@]} -gt 0 ]; then
            MERGED_VCF="${merged_candidates[0]}"
            break
        fi
    done
fi

if [ -n "${MERGED_VCF}" ]; then
    LOCAL_MERGED_VCF="combined_for_pca.vcf.gz"
    if [ "${MERGED_VCF}" != "${PWD}/${LOCAL_MERGED_VCF}" ]; then
        if ln -sf "${MERGED_VCF}" "${LOCAL_MERGED_VCF}"; then
            if [ -f "${MERGED_VCF}.tbi" ]; then
                ln -sf "${MERGED_VCF}.tbi" "${LOCAL_MERGED_VCF}.tbi"
            fi
            COMBINED_VCF="${LOCAL_MERGED_VCF}"
        else
            log_warn "Failed to link merged VCF; using source path instead."
            COMBINED_VCF="${MERGED_VCF}"
        fi
    else
        COMBINED_VCF="${MERGED_VCF}"
    fi
    log_info "Using merged VCF for PCA: ${COMBINED_VCF}"
    if [[ "${COMBINED_VCF}" == *.vcf.gz ]] && [ ! -f "${COMBINED_VCF}.tbi" ]; then
        log_info "Indexing merged VCF: ${COMBINED_VCF}"
        "${BCFTOOLS_BIN}" index -t "${COMBINED_VCF}"
    fi
else
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
fi

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

DUPLICATE_SAMPLES_FILE=""
if [ "${DUPLICATE_MODE}" != "off" ]; then
    log_info "Estimating KING kinship for duplicate detection (threshold ${DUPLICATE_THRESHOLD})"
    KING_PREFIX="king_duplicates"
    "${PLINK2_BIN}" \
        --pfile qc \
        --king \
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
    --pca 10 var-wts \
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
