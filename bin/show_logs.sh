#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_ROOT="${PIPELINE_ROOT:-$(cd "${SCRIPT_DIR}/.." && pwd)}"

usage() {
    echo "Usage: ${0##*/} <dataset>"
    exit 1
}

[ $# -ge 1 ] || usage
DATASET="$1"

# Source config and logging helpers
if [ -f "${PIPELINE_ROOT}/config/pipeline_config.sh" ]; then
    # shellcheck source=../config/pipeline_config.sh
    source "${PIPELINE_ROOT}/config/pipeline_config.sh"
fi
if [ -f "${PIPELINE_ROOT}/lib/logging.sh" ]; then
    # shellcheck source=../lib/logging.sh
    source "${PIPELINE_ROOT}/lib/logging.sh"
fi

# Resolve key paths
log_base="$(resolve_log_root "${DATASET}" "")"
pipeline_dir="$(resolve_log_root "${DATASET}" "pipeline")"
slurm_dir="$(resolve_log_root "${DATASET}" "slurm")"
artifacts_dir="$(resolve_log_root "${DATASET}" "artifacts")"
state_dir="${MASTER_LOG_DIR:-${LOG_BASE_PATH%/}/${DATASET}}"

cat <<EOF
Dataset: ${DATASET}
Log root: ${log_base}
- pipeline: ${pipeline_dir}
- slurm:    ${slurm_dir}
- artifacts:${artifacts_dir}
- state:    ${state_dir}

Recent pipeline logs:
$(ls -1t "${pipeline_dir}" 2>/dev/null | head -10 || true)

Recent slurm logs:
$(ls -1t "${slurm_dir}" 2>/dev/null | head -10 || true)
EOF
