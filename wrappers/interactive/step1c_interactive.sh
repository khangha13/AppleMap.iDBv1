#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [ -n "${PIPELINE_ROOT:-}" ]; then
    if [ -d "${PIPELINE_ROOT}" ]; then
        PIPELINE_ROOT="$(cd "${PIPELINE_ROOT}" && pwd)"
    else
        echo "[step1c_interactive] ⚠️  Provided PIPELINE_ROOT (${PIPELINE_ROOT}) does not exist; falling back to wrapper-relative path." >&2
        PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
    fi
else
    PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
fi
export PIPELINE_ROOT
MODULE_DIR="${PIPELINE_ROOT}/modules/step1c"

source "${PIPELINE_ROOT}/lib/logging.sh"
source "${PIPELINE_ROOT}/lib/validation.sh"
source "${PIPELINE_ROOT}/lib/interactive_common.sh"
source "${PIPELINE_ROOT}/config/pipeline_config.sh"
source "${MODULE_DIR}/lib/functions.sh"

init_logging "step1c_interactive" "interactive"

main() {
    echo "\n🧬 STEP 1C INTERACTIVE LAUNCHER"
    echo "==============================="

    local dataset_name
    dataset_name=$(prompt_dataset_name)

    local default_rdm="$(get_rdm_datasets_path)/${dataset_name}"
    local rdm_base_path
    rdm_base_path=$(prompt_directory "RDM base path" "${default_rdm}")

    if [ ! -d "${rdm_base_path}/7.Consolidated_VCF" ]; then
        echo "❌ Required directory missing: ${rdm_base_path}/7.Consolidated_VCF"
        exit 1
    fi

    validate_step1c_inputs "${rdm_base_path}" || true

    echo "\n📂 Dataset: ${dataset_name}"
    echo "📁 RDM path: ${rdm_base_path}"

    local default_output
    default_output=$(default_step1c_output_dir "${rdm_base_path}")
    read -p "Output directory [${default_output}]: " STEP1C_OUTPUT_DIR
    STEP1C_OUTPUT_DIR=${STEP1C_OUTPUT_DIR:-${default_output}}
    mkdir -p "${STEP1C_OUTPUT_DIR}" || true

    local gene_map
    read -p "Gene map file (optional): " gene_map
    if [ -n "${gene_map}" ] && [ ! -f "${gene_map}" ]; then
        echo "❌ Gene map file not found: ${gene_map}"
        exit 1
    fi

    echo "\n🗂️  Joint VCFs discovered:"
    find "${rdm_base_path}/7.Consolidated_VCF" -maxdepth 1 -name '*.vcf.gz' -type f | head -n 10 | sed 's/^/  • /'

    if confirm_action "Proceed with Step 1C submission?"; then
        GENE_MAP_FILE="${gene_map}" STEP1C_OUTPUT_DIR="${STEP1C_OUTPUT_DIR}" \
            bash "${PIPELINE_ROOT}/wrappers/sbatch/step1c_submit.sh" "${dataset_name}" "${rdm_base_path}"
    else
        echo "Operation cancelled."
    fi
}

main "$@"
