#!/bin/bash
# =============================================================================
# Optional environment bootstrap for resolving PIPELINE_ROOT on HPC clusters
# =============================================================================

# Allow admins/users to override the file that stores persistent exports.
if [ -z "${PIPELINE_ENV_FILE:-}" ]; then
    PIPELINE_ENV_FILE="${HOME%/}/.gatk_pipeline_env"
fi

if [ -f "${PIPELINE_ENV_FILE}" ]; then
    # shellcheck disable=SC1090
    source "${PIPELINE_ENV_FILE}"
fi

# Provide sane defaults that other scripts can reference.
PIPELINE_DIR_NAME="${PIPELINE_DIR_NAME:-GATK_Pipeline_KH_v1}"

if [ -z "${PIPELINE_HOME_CANDIDATE:-}" ] && [ -n "${HOME:-}" ]; then
    PIPELINE_HOME_CANDIDATE="${HOME%/}/${PIPELINE_DIR_NAME}"
fi

export PIPELINE_ENV_FILE PIPELINE_DIR_NAME PIPELINE_HOME_CANDIDATE

