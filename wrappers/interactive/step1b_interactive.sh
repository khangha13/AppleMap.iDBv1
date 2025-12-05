#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [ -n "${PIPELINE_ROOT:-}" ]; then
    if [ -d "${PIPELINE_ROOT}" ]; then
        PIPELINE_ROOT="$(cd "${PIPELINE_ROOT}" && pwd)"
    else
        echo "[step1b_interactive] ⚠️  Provided PIPELINE_ROOT (${PIPELINE_ROOT}) does not exist; falling back to wrapper-relative path." >&2
        PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
    fi
else
    PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
fi
export PIPELINE_ROOT
MODULE_DIR="${PIPELINE_ROOT}/modules/step1b"

source "${PIPELINE_ROOT}/lib/logging.sh"
source "${PIPELINE_ROOT}/lib/validation.sh"
source "${PIPELINE_ROOT}/lib/interactive_common.sh"
source "${PIPELINE_ROOT}/config/pipeline_config.sh"
source "${MODULE_DIR}/lib/functions.sh"

init_logging "step1b_interactive" "interactive"

main() {
    echo "\n🧬 STEP 1B INTERACTIVE LAUNCHER"
    echo "==============================="

    local dataset_name
    dataset_name=$(prompt_dataset_name)

    local default_rdm="$(get_rdm_datasets_path)/${dataset_name}"
    local rdm_base_path
    rdm_base_path=$(prompt_directory "RDM base path" "${default_rdm}")

    if [ ! -d "${rdm_base_path}/5.Individual_VCF" ]; then
        echo "❌ Required directory missing: ${rdm_base_path}/5.Individual_VCF"
        exit 1
    fi

    echo "\n📂 Dataset: ${dataset_name}"
    echo "📁 RDM path: ${rdm_base_path}"

    local reference_genome
    reference_genome="$(get_reference_fasta)"
    if [ ! -f "${reference_genome}" ]; then
        echo "❌ Reference genome not found: ${reference_genome}"
        exit 1
    fi

    echo "\n🧬 Chromosomes detected:"
    get_chromosome_list "${reference_genome}" | sed 's/^/  • /'

    if confirm_action "Proceed with Step 1B submission?"; then
        bash "${PIPELINE_ROOT}/wrappers/sbatch/step1b_submit.sh" "${dataset_name}" "${rdm_base_path}"
    else
        echo "Operation cancelled."
    fi
}

main "$@"
