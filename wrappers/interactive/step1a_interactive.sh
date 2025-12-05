#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [ -n "${PIPELINE_ROOT:-}" ]; then
    if [ -d "${PIPELINE_ROOT}" ]; then
        PIPELINE_ROOT="$(cd "${PIPELINE_ROOT}" && pwd)"
    else
        echo "[step1a_interactive] ⚠️  Provided PIPELINE_ROOT (${PIPELINE_ROOT}) does not exist; falling back to wrapper-relative path." >&2
        PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
    fi
else
    PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
fi
export PIPELINE_ROOT

source "${PIPELINE_ROOT}/lib/logging.sh"
source "${PIPELINE_ROOT}/lib/validation.sh"
source "${PIPELINE_ROOT}/lib/interactive_common.sh"
source "${PIPELINE_ROOT}/config/pipeline_config.sh"

init_logging "step1a_interactive" "interactive"

main() {
    echo "\n🧬 STEP 1A INTERACTIVE LAUNCHER"
    echo "==============================="

    # Optional flags:
    #   --sample=NAME        Run a single sample (basename of FASTQs)
    #   --sample-list=PATH   Run samples from a custom list (one per line)
    local sample_override=""
    local sample_list_override=""
    for arg in "$@"; do
        case "${arg}" in
            --sample=*)
                sample_override="${arg#--sample=}"
                ;;
            --sample-list=*)
                sample_list_override="${arg#--sample-list=}"
                ;;
            *)
                ;;
        esac
    done

    local dataset_name
    dataset_name=$(prompt_dataset_name)

    local default_rdm="$(get_rdm_datasets_path)/${dataset_name}"
    local rdm_base_path
    rdm_base_path=$(prompt_directory "RDM base path (If correct, press 'enter' to confirm)" "${default_rdm}")

    if [ ! -d "${rdm_base_path}/1.FASTQ" ]; then
        echo "❌ Expected directory missing: ${rdm_base_path}/1.FASTQ"
        exit 1
    fi

    echo "\n📂 Dataset: ${dataset_name}"
    echo "📁 RDM path: ${rdm_base_path}"
    echo "ℹ️  When prompted next, press 'y' to submit or 'n' to cancel."
    if [ -n "${sample_override}" ]; then
        echo "🎯 Single-sample mode: ${sample_override}"
    fi
    if [ -n "${sample_list_override}" ]; then
        echo "📝 Using sample list: ${sample_list_override}"
    fi

    validate_sample_files "${rdm_base_path}/1.FASTQ" || true

    if confirm_action "Proceed with Step 1A submission?"; then
        extra_args=()
        if [ -n "${sample_override}" ]; then
            extra_args+=(--sample "${sample_override}")
        fi
        if [ -n "${sample_list_override}" ]; then
            extra_args+=(--sample-list "${sample_list_override}")
        fi
        bash "${PIPELINE_ROOT}/wrappers/sbatch/step1a_submit.sh" "${dataset_name}" "${rdm_base_path}" "${extra_args[@]}"
    else
        echo "Operation cancelled."
    fi
}

main "$@"
