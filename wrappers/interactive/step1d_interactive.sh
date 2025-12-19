#!/bin/bash
# =============================================================================
# STEP 1D INTERACTIVE WRAPPER
# =============================================================================
# Prompts the user for a working directory (containing Chr00–Chr17 VCFs),
# confirms the selection, and launches master_vcf_analysis.sh with the
# appropriate environment variables set.
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [ -n "${PIPELINE_ROOT:-}" ]; then
    if [ -d "${PIPELINE_ROOT}" ]; then
        PIPELINE_ROOT="$(cd "${PIPELINE_ROOT}" && pwd)"
    else
        echo "[step1d_interactive] ⚠️  Provided PIPELINE_ROOT (${PIPELINE_ROOT}) does not exist; falling back to wrapper-relative path." >&2
        PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
    fi
else
    PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
fi
export PIPELINE_ROOT
MODULE_DIR="${PIPELINE_ROOT}/modules/step1d"
if [ -f "${PIPELINE_ROOT}/config/pipeline_config.sh" ]; then
    source "${PIPELINE_ROOT}/config/pipeline_config.sh"
fi

usage() {
    cat <<'EOF'
Usage: step1d_interactive.sh [--dir=PATH] [--vcf=NAME[,NAME...]] [--beagle] [--pca-only] [--remove-relatives] [--dry-run] [--help]

Options:
  --dir=PATH      Full path to directory containing VCF files (skips prompt).
  --vcf=LIST      Comma-separated list of VCF prefixes (e.g. Chr00,chr01).
                  Matching is case-insensitive; extensions default to .vcf.gz.
  --beagle        Pass --beagle to master_vcf_analysis.sh (Beagle-imputed VCFs).
  --qc            Run QC-only mode (metrics/plots, default).
  --PCA           Run PCA-only mode (expects merged or per-chrom VCFs).
  --duplicate-check
                  Run KING duplicate detection only (writes tables).
  --remove-relatives
                  (Requires --PCA) Remove close relatives before PCA using KING 0.125.
  --dry-run       Pass --dry-run to master_vcf_analysis.sh (preview without creating files).
  --help          Show this message and exit.
EOF
}

VCF_OVERRIDE=""
BEAGLE_FLAG=false
DRY_RUN_FLAG=false
DIR_OVERRIDE=""
if [ -n "${STEP1D_EXPECTED_CHROMS:-}" ]; then
    EXPECTED_CHROM_COUNT="${STEP1D_EXPECTED_CHROMS}"
else
    EXPECTED_CHROM_COUNT=18
fi
EXPECTED_CHROM_DESC="Chr00–Chr17"
MODE="qc"
MODE_SET=false
REMOVE_REL_FLAG=false

declare -a AVAILABLE_PRIMARY_VCFS=()
declare -a AVAILABLE_OTHER_VCFS=()
declare -a FILTERED_VCF_CANDIDATES=()
declare -a AUTO_VCF_BASENAMES=()
declare -a SELECTED_VCF_BASENAMES=()
declare -a MISSING_VCF_REQUESTS=()
declare -a EFFECTIVE_VCF_BASENAMES=()
FALLBACK_MODE=false
USING_FILTERED_VCFS=false
INFERRED_VCF_PATTERN=""

for arg in "$@"; do
case "${arg}" in
        --dir=*)
            DIR_OVERRIDE="${arg#--dir=}"
            ;;
        --vcf=*)
            VCF_OVERRIDE="${arg#--vcf=}"
            ;;
        --beagle)
            BEAGLE_FLAG=true
            ;;
        --qc)
            if ${MODE_SET}; then
                echo "❌ Multiple modes provided; choose one of --qc, --PCA, or --duplicate-check." >&2
                exit 1
            fi
            MODE="qc"
            MODE_SET=true
            ;;
        --PCA|--pca|--pca-only)
            if ${MODE_SET}; then
                echo "❌ Multiple modes provided; choose one of --qc, --PCA, or --duplicate-check." >&2
                exit 1
            fi
            MODE="pca"
            MODE_SET=true
            ;;
        --duplicate-check)
            if ${MODE_SET}; then
                echo "❌ Multiple modes provided; choose one of --qc, --PCA, or --duplicate-check." >&2
                exit 1
            fi
            MODE="duplicate-check"
            MODE_SET=true
            ;;
        --remove-relatives)
            REMOVE_REL_FLAG=true
            ;;
        --dry-run|-n)
            DRY_RUN_FLAG=true
            ;;
        --help|-h)
            usage
            exit 0
            ;;
        *)
            echo "❌ Unknown option: ${arg}" >&2
            usage
            exit 1
            ;;
    esac
done

if [ "${BEAGLE_FLAG}" = true ] && [ -z "${STEP1D_EXPECTED_CHROMS:-}" ]; then
    EXPECTED_CHROM_COUNT=17
fi
if [ "${BEAGLE_FLAG}" = true ] && [ "${EXPECTED_CHROM_COUNT}" -eq 17 ]; then
    EXPECTED_CHROM_DESC="Chr01–Chr17 (Beagle mode)"
elif [ -n "${STEP1D_EXPECTED_CHROMS:-}" ]; then
    EXPECTED_CHROM_DESC="custom configuration"
fi

if [ "${REMOVE_REL_FLAG}" = true ] && [ "${MODE}" != "pca" ]; then
    echo "❌ --remove-relatives requires --PCA." >&2
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MASTER_SCRIPT="${MODULE_DIR}/templates/master_vcf_analysis.sh"

if [ ! -f "${MASTER_SCRIPT}" ]; then
    echo "❌ Cannot locate master_vcf_analysis.sh at ${MASTER_SCRIPT}" >&2
    exit 1
fi

prompt_for_directory() {
    local input_path=""
    while true; do
        read -r -p "Enter the full path to the working directory: " input_path
        if [ -z "${input_path}" ]; then
            echo "⚠️  Path cannot be empty. Please try again."
            continue
        fi
        if [ ! -d "${input_path}" ]; then
            echo "❌ Directory does not exist: ${input_path}"
            continue
        fi
        break
    done
    echo "${input_path}"
}

confirm_selection() {
    local prompt="$1"
    local response=""
    while true; do
        read -r -p "${prompt} (y/n): " response
        case "${response}" in
            [Yy]*) return 0 ;;
            [Nn]*) return 1 ;;
            *) echo "Please answer y or n." ;;
        esac
    done
}

summarise_directory() {
    local dir_path="$1"
    echo ""
    echo "📂 Selected working directory: ${dir_path}"
    AVAILABLE_PRIMARY_VCFS=()
    AVAILABLE_OTHER_VCFS=()
    FILTERED_VCF_CANDIDATES=()
    USING_FILTERED_VCFS=false

    while IFS= read -r -d '' path; do
        local base
        base=$(basename "${path}")
        case "${base}" in
            *_snps.vcf.gz|*_snps.vcf|*_snps.vcf.fz|*_snp.vcf.gz|*_snp.vcf.fz)
                FILTERED_VCF_CANDIDATES+=("${base}")
                continue
            ;;
        esac
        if [[ "${base}" =~ ^Chr(0[0-9]|1[0-7])\.vcf\.gz$ ]]; then
            AVAILABLE_PRIMARY_VCFS+=("${base}")
        else
            AVAILABLE_OTHER_VCFS+=("${base}")
        fi
    done < <(find "${dir_path}" -maxdepth 1 -type f \( -name "*.vcf" -o -name "*.vcf.gz" \) -print0)

    if [ ${#AVAILABLE_PRIMARY_VCFS[@]} -gt 0 ]; then
        mapfile -t AVAILABLE_PRIMARY_VCFS < <(printf '%s\n' "${AVAILABLE_PRIMARY_VCFS[@]}" | sort)
    fi
    if [ ${#AVAILABLE_OTHER_VCFS[@]} -gt 0 ]; then
        mapfile -t AVAILABLE_OTHER_VCFS < <(printf '%s\n' "${AVAILABLE_OTHER_VCFS[@]}" | sort)
    fi
    if [ ${#FILTERED_VCF_CANDIDATES[@]} -gt 0 ]; then
        mapfile -t FILTERED_VCF_CANDIDATES < <(printf '%s\n' "${FILTERED_VCF_CANDIDATES[@]}" | sort)
    fi

    local primary_count=${#AVAILABLE_PRIMARY_VCFS[@]}
    local other_count=${#AVAILABLE_OTHER_VCFS[@]}
    local filtered_count=${#FILTERED_VCF_CANDIDATES[@]}

    if [ "${primary_count}" -eq 0 ] && [ "${filtered_count}" -gt 0 ]; then
        AVAILABLE_PRIMARY_VCFS=("${FILTERED_VCF_CANDIDATES[@]}")
        primary_count=${#AVAILABLE_PRIMARY_VCFS[@]}
        USING_FILTERED_VCFS=true
    fi

    if [ "${primary_count}" -gt 0 ]; then
        if [ "${USING_FILTERED_VCFS}" = true ]; then
            echo "✅ Found ${primary_count} filtered per-chromosome VCF(s) ending with *_snps.vcf.gz."
        else
            echo "✅ Found ${primary_count} per-chromosome VCF(s) matching Chr??.vcf.gz."
        fi
    else
        echo "⚠️  No files matching Chr??.vcf.gz were found (excluding *_snps.vcf.gz)."
    fi

    if [ "${other_count}" -gt 0 ]; then
        echo "ℹ️  Additional VCF candidates detected:"
        for base in "${AVAILABLE_OTHER_VCFS[@]}"; do
            echo "    • ${base}"
        done
    else
        echo "ℹ️  No other VCF files detected."
    fi

    if [ "${filtered_count}" -gt 0 ] && [ "${USING_FILTERED_VCFS}" = false ]; then
        echo "ℹ️  Filtered *_snps VCF files detected (${filtered_count}); they will be used automatically if raw Chr?? files are absent."
    elif [ "${USING_FILTERED_VCFS}" = true ]; then
        echo "ℹ️  Using filtered *_snps VCF files (${primary_count}) because no raw Chr??.vcf(.gz) files were found."
    fi

    echo "ℹ️  This pipeline expects ${EXPECTED_CHROM_COUNT} per-chromosome VCFs (${EXPECTED_CHROM_DESC}). Detected ${primary_count}."
    if [ "${primary_count}" -ne "${EXPECTED_CHROM_COUNT}" ]; then
        echo "   Tip: export VCF_PATTERN to match your naming scheme if needed."
    fi
}

trim() {
    local value="$*"
    value="${value#"${value%%[![:space:]]*}"}"
    value="${value%"${value##*[![:space:]]}"}"
    printf '%s' "${value}"
}

validate_requested_vcfs() {
    local raw_list="$1"
    IFS=',' read -r -a requested <<< "${raw_list}"
    SELECTED_VCF_BASENAMES=()
    MISSING_VCF_REQUESTS=()
    for item in "${requested[@]}"; do
        local trimmed
        trimmed="$(trim "${item}")"
        [ -z "${trimmed}" ] && continue
        local candidate="${trimmed}"
        if [[ "${candidate}" =~ \.[Vv][Cc][Ff](\.gz)?$ ]]; then
            :
        else
            candidate="${candidate}.vcf.gz"
        fi
        local match
        match=$(find "${WORK_DIR_ABS}" -maxdepth 1 -type f -iname "${candidate}" 2>/dev/null | head -n 1 || true)
        if [ -z "${match}" ] && [[ ! "${trimmed}" =~ \.[Vv][Cc][Ff](\.gz)?$ ]]; then
            local fallback="${trimmed}.vcf"
            match=$(find "${WORK_DIR_ABS}" -maxdepth 1 -type f -iname "${fallback}" 2>/dev/null | head -n 1 || true)
        fi
        if [ -n "${match}" ]; then
            local base
            base=$(basename "${match}")
            local exists=false
            for selected in "${SELECTED_VCF_BASENAMES[@]}"; do
                if [ "${selected}" = "${base}" ]; then
                    exists=true
                    break
                fi
            done
            if [ "${exists}" = false ]; then
                SELECTED_VCF_BASENAMES+=("${base}")
            fi
        else
            MISSING_VCF_REQUESTS+=("${trimmed}")
        fi
    done
}

infer_vcf_pattern() {
    local -a filenames=("$@")
    if [ ${#filenames[@]} -eq 0 ]; then
        return 1
    fi

    local base
    local prefix=""
    local chr_literal=""
    local suffix=""

    for base in "${filenames[@]}"; do
        if [[ "${base}" =~ ^(.*?)([Cc]hr)([0-9]{2})(.*)$ ]]; then
            prefix="${BASH_REMATCH[1]}"
            chr_literal="${BASH_REMATCH[2]}"
            suffix="${BASH_REMATCH[4]}"
            break
        fi
    done

    if [ -z "${chr_literal}" ]; then
        return 1
    fi

    for base in "${filenames[@]}"; do
        if [[ "${base}" =~ ^(.*?)([Cc]hr)([0-9]{2})(.*)$ ]]; then
            if [ "${BASH_REMATCH[1]}" != "${prefix}" ] || [ "${BASH_REMATCH[2]}" != "${chr_literal}" ] || [ "${BASH_REMATCH[4]}" != "${suffix}" ]; then
                return 1
            fi
        else
            return 1
        fi
    done

    printf '%s%s%%02d%s\n' "${prefix}" "${chr_literal}" "${suffix}"
}

echo "🧬 STEP 1D INTERACTIVE LAUNCHER"
echo "==============================="
echo ""
echo "This helper will collect the required path information and run:"
echo "  ${MASTER_SCRIPT##*/}"
echo ""

if [ -n "${DIR_OVERRIDE}" ]; then
    # Use provided directory
    WORK_DIR_INPUT="${DIR_OVERRIDE}"
    if [ ! -d "${WORK_DIR_INPUT}" ]; then
        echo "❌ Directory does not exist: ${WORK_DIR_INPUT}" >&2
        exit 1
    fi
else
    # Prompt for directory
    WORK_DIR_INPUT="$(prompt_for_directory)"
fi

# Resolve to absolute path without relying on GNU readlink -f
if ! WORK_DIR_ABS="$(cd "${WORK_DIR_INPUT}" && pwd)"; then
    echo "❌ Unable to resolve absolute path for ${WORK_DIR_INPUT}" >&2
    exit 1
fi

summarise_directory "${WORK_DIR_ABS}"
echo ""

if [ -z "${VCF_OVERRIDE}" ] && [ ${#AVAILABLE_PRIMARY_VCFS[@]} -eq 0 ] && [ ${#AVAILABLE_OTHER_VCFS[@]} -gt 0 ]; then
    FALLBACK_MODE=true
    AUTO_VCF_BASENAMES=("${AVAILABLE_OTHER_VCFS[@]}")
else
    FALLBACK_MODE=false
    AUTO_VCF_BASENAMES=()
fi

if [ -n "${VCF_OVERRIDE}" ]; then
    echo "🔍 Validating requested VCF subset..."
    validate_requested_vcfs "${VCF_OVERRIDE}"
    echo ""
    echo "🎯 Requested VCF files:"
    for base in "${SELECTED_VCF_BASENAMES[@]}"; do
        echo "   ✅ ${base}"
    done
    if [ ${#MISSING_VCF_REQUESTS[@]} -gt 0 ]; then
        for missing in "${MISSING_VCF_REQUESTS[@]}"; do
            echo "   ❌ ${missing} (not found)"
        done
        echo "⚠️  Missing entries will be skipped."
    fi
    if [ ${#SELECTED_VCF_BASENAMES[@]} -eq 0 ]; then
        echo "❌ None of the requested VCF files were found in ${WORK_DIR_ABS}. Aborting."
        exit 1
    fi
fi

echo ""
EFFECTIVE_VCF_BASENAMES=()
if [ ${#SELECTED_VCF_BASENAMES[@]} -gt 0 ]; then
    EFFECTIVE_VCF_BASENAMES=("${SELECTED_VCF_BASENAMES[@]}")
    confirm_message="Proceed with Step 1D for ${#EFFECTIVE_VCF_BASENAMES[@]} requested VCF file(s)?"
elif [ "${FALLBACK_MODE}" = "true" ]; then
    EFFECTIVE_VCF_BASENAMES=("${AUTO_VCF_BASENAMES[@]}")
    echo "⚠️  No Chr00–Chr17 matches were found. The following VCF files will be used if you continue:"
    for base in "${EFFECTIVE_VCF_BASENAMES[@]}"; do
        echo "   • ${base}"
    done
    confirm_message="Proceed with Step 1D using these ${#EFFECTIVE_VCF_BASENAMES[@]} VCF file(s)?"
else
    EFFECTIVE_VCF_BASENAMES=("${AVAILABLE_PRIMARY_VCFS[@]}")
    detected_count=${#EFFECTIVE_VCF_BASENAMES[@]}
    if [ "${detected_count}" -eq 0 ]; then
        confirm_message="No ${EXPECTED_CHROM_DESC} VCFs detected. Continue anyway?"
    elif [ "${detected_count}" -ne "${EXPECTED_CHROM_COUNT}" ]; then
        confirm_message="Detected ${detected_count} of ${EXPECTED_CHROM_COUNT} expected VCFs (${EXPECTED_CHROM_DESC}). Continue anyway?"
    else
        confirm_message="Proceed with Step 1D using the detected ${EXPECTED_CHROM_DESC} VCF files?"
    fi
fi

if ! confirm_selection "${confirm_message}"; then
    echo "Operation cancelled."
    exit 0
fi

PATTERN_NOTE=""
INFERRED_VCF_PATTERN=""
if [ -z "${VCF_PATTERN:-}" ] && [ ${#EFFECTIVE_VCF_BASENAMES[@]} -gt 0 ]; then
    if INFERRED_VCF_PATTERN="$(infer_vcf_pattern "${EFFECTIVE_VCF_BASENAMES[@]}")"; then
        export VCF_PATTERN="${INFERRED_VCF_PATTERN}"
    else
        INFERRED_VCF_PATTERN=""
    fi
fi

if [ -n "${VCF_PATTERN:-}" ]; then
    if [ -n "${INFERRED_VCF_PATTERN}" ] && [ "${VCF_PATTERN}" = "${INFERRED_VCF_PATTERN}" ]; then
        PATTERN_NOTE="    Pattern:  ${VCF_PATTERN} (auto-inferred)"
    else
        PATTERN_NOTE="    Pattern:  ${VCF_PATTERN} (pre-set)"
    fi
else
    PATTERN_NOTE="    Pattern:  (default Chr%02d.vcf.gz; set VCF_PATTERN manually to override)"
fi

echo ""
echo "🚀 Launching master_vcf_analysis.sh..."
echo "    VCF_DIR=${WORK_DIR_ABS}"
echo "    WORK_DIR=${WORK_DIR_ABS}"
echo "${PATTERN_NOTE}"
if [ "${BEAGLE_FLAG}" = "true" ]; then
    echo "    Mode:    --beagle (Beagle-imputed inputs)"
fi
if [ "${DRY_RUN_FLAG}" = "true" ]; then
    echo "    Mode:    --dry-run (preview mode, no files created)"
fi
case "${MODE}" in
    qc)
        echo "    Mode:    QC (metrics + plots + AF distributions)"
        ;;
    pca)
        echo "    Mode:    PCA-only (merged VCF preferred)"
        if [ "${REMOVE_REL_FLAG}" = "true" ]; then
            echo "             └─ Removing close relatives (KING cutoff 0.125)"
        fi
        ;;
    duplicate-check)
        echo "    Mode:    Duplicate-check only (KING duplicates table)"
        ;;
esac
echo ""

if [ -z "${VCF_PATTERN:-}" ] && [ -z "${INFERRED_VCF_PATTERN}" ]; then
    echo "⚠️  Could not infer a VCF filename pattern automatically. If your files are not named like Chr%02d.vcf.gz, export VCF_PATTERN before rerunning."
fi

export VCF_DIR="${WORK_DIR_ABS}"
export WORK_DIR="${WORK_DIR_ABS}"
export R_SCRIPTS_DIR="${MODULE_DIR}/Rscripts"
if [ ${#EFFECTIVE_VCF_BASENAMES[@]} -gt 0 ]; then
    export VCF_INCLUDE_FILENAMES="${EFFECTIVE_VCF_BASENAMES[*]}"
else
    unset VCF_INCLUDE_FILENAMES
fi

if [ -z "${TMPDIR:-}" ]; then
    TMPDIR="$(mktemp -d "${SCRATCH_BASE_PATH%/}/step1d_interactive.XXXXXX")"
    export TMPDIR
    trap 'rm -rf "${TMPDIR}"' EXIT
fi

declare -a MASTER_ARGS=()
if [ "${BEAGLE_FLAG}" = "true" ]; then
    MASTER_ARGS+=("--beagle")
fi
if [ "${DRY_RUN_FLAG}" = "true" ]; then
    MASTER_ARGS+=("--dry-run")
fi
case "${MODE}" in
    qc) MASTER_ARGS+=("--qc") ;;
    pca) MASTER_ARGS+=("--PCA") ;;
    duplicate-check) MASTER_ARGS+=("--duplicate-check") ;;
esac
if [ "${REMOVE_REL_FLAG}" = "true" ]; then
    MASTER_ARGS+=("--remove-relatives")
fi

if ! bash "${MASTER_SCRIPT}" "${MASTER_ARGS[@]}"; then
    echo "❌ master_vcf_analysis.sh exited with an error." >&2
    exit 1
fi

echo ""
echo "🎉 Interactive Step 1D complete."
