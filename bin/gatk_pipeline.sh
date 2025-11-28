#!/bin/bash -l
#SBATCH --job-name=GATK_master_script
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=336:00:00
#SBATCH --account=a_qaafi_cas
#SBATCH --partition=general
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

BOOTSTRAP_PATH="${SCRIPT_DIR}/../lib/env_bootstrap.sh"
if [ -f "${BOOTSTRAP_PATH}" ]; then
    # shellcheck source=../lib/env_bootstrap.sh
    source "${BOOTSTRAP_PATH}"
fi
unset BOOTSTRAP_PATH

resolve_pipeline_root() {
    local candidate resolved
    local -a search_paths=()

if [ -n "${PIPELINE_ROOT:-}" ]; then
        search_paths+=( "${PIPELINE_ROOT}" )
    fi
    if [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
        search_paths+=( "${SLURM_SUBMIT_DIR}" "${SLURM_SUBMIT_DIR}/.." )
    fi
    search_paths+=( "${SCRIPT_DIR}/.." "${SCRIPT_DIR}" )

    local home_fallback=""
    if [ -n "${PIPELINE_HOME_CANDIDATE:-}" ]; then
        home_fallback="${PIPELINE_HOME_CANDIDATE}"
    elif [ -n "${HOME:-}" ]; then
        home_fallback="${HOME%/}/${PIPELINE_DIR_NAME:-GATK_Pipeline_KH_v1}"
    fi

    if [ -n "${home_fallback}" ]; then
        search_paths+=( "${home_fallback}" )
    fi

    for candidate in "${search_paths[@]}"; do
        [ -n "${candidate}" ] || continue
        if resolved="$(cd "${candidate}" 2>/dev/null && pwd)"; then
            if [ -f "${resolved}/lib/logging.sh" ]; then
                PIPELINE_ROOT="${resolved}"
export PIPELINE_ROOT
                return 0
            fi
        fi
    done

    echo "[gatk_pipeline.sh] ❌ Unable to locate pipeline root. Please rerun from the repository root or set PIPELINE_ROOT explicitly." >&2
    exit 1
}

resolve_pipeline_root

source "${PIPELINE_ROOT}/lib/logging.sh"
source "${PIPELINE_ROOT}/lib/validation.sh"
source "${PIPELINE_ROOT}/lib/slurm.sh"
source "${PIPELINE_ROOT}/lib/pipeline_common.sh"
source "${PIPELINE_ROOT}/config/pipeline_config.sh"

MASTER_RUN_STAMP="${MASTER_RUN_STAMP:-$(date +%Y%m%d_%H%M%S)}"
MASTER_LOG_DIR="${MASTER_LOG_DIR:-}"
MASTER_LOG_FILE="${MASTER_LOG_FILE:-}"
MASTER_LOG_DATASET="${MASTER_LOG_DATASET:-}"
MASTER_LOG_STDOUT="${MASTER_LOG_STDOUT:-}"
MASTER_LOG_STDERR="${MASTER_LOG_STDERR:-}"
MASTER_LOG_INITIALIZED="${MASTER_LOG_INITIALIZED:-}"
MASTER_IO_REDIRECTED="${MASTER_IO_REDIRECTED:-}"

ensure_master_logging() {
    local dataset_name="$1"
    local mode="${2:-full}"

    if [ -z "${dataset_name:-}" ]; then
        return 0
    fi

    if [ -z "${MASTER_RUN_STAMP:-}" ]; then
        MASTER_RUN_STAMP="$(date +%Y%m%d_%H%M%S)"
        export MASTER_RUN_STAMP
    fi

    if [ "${MASTER_LOG_DATASET:-}" != "${dataset_name}" ] || [ -z "${MASTER_LOG_DIR:-}" ]; then
        local log_base="${LOG_BASE_PATH:-${SCRATCH_BASE_PATH%/}/logs}"
        MASTER_LOG_DIR="${log_base%/}/${dataset_name}_${MASTER_RUN_STAMP}"
        MASTER_LOG_DATASET="${dataset_name}"
    fi

    mkdir -p "${MASTER_LOG_DIR}"
    local log_suffix="${SLURM_JOB_ID:-${MASTER_RUN_STAMP}}"
    MASTER_LOG_STDOUT="${MASTER_LOG_DIR}/GATK_master_script_${log_suffix}.output"
    MASTER_LOG_STDERR="${MASTER_LOG_DIR}/GATK_master_script_${log_suffix}.error"
    export MASTER_LOG_DIR MASTER_LOG_STDOUT MASTER_LOG_STDERR MASTER_LOG_DATASET

    if [ "${mode}" = "prepare" ]; then
        return 0
    fi

    if [ -z "${MASTER_LOG_INITIALIZED:-}" ]; then
        local previous_override="${PIPELINE_LOG_DIR_OVERRIDE:-}"
        PIPELINE_LOG_DIR_OVERRIDE="${MASTER_LOG_DIR}"
        init_logging "gatk_master" "pipeline"
        MASTER_LOG_FILE="${LOG_FILE}"
        MASTER_LOG_INITIALIZED="true"
        export MASTER_LOG_FILE MASTER_LOG_INITIALIZED
        if [ -n "${previous_override}" ]; then
            PIPELINE_LOG_DIR_OVERRIDE="${previous_override}"
        else
            unset PIPELINE_LOG_DIR_OVERRIDE
        fi
        log_info "Master logging directory: ${MASTER_LOG_DIR}"
    fi

    if [ -z "${MASTER_IO_REDIRECTED:-}" ]; then
        exec > >(tee -a "${MASTER_LOG_STDOUT}") 2> >(tee -a "${MASTER_LOG_STDERR}" >&2)
        MASTER_IO_REDIRECTED="true"
        export MASTER_IO_REDIRECTED
    fi
}

#
# Master loop utilities (used when submitted via sbatch or when --step=auto with dataset)
#
poll_interval_seconds="${PIPELINE_MASTER_POLL_SECS:-${MONITOR_INTERVAL:-60}}"
master_stall_limit="${PIPELINE_MASTER_STALL_POLLS:-10}"

cancel_master_job() {
    local reason="${1:-Unspecified reason}"
    log_error "Cancelling master job: ${reason}"
    if [ -n "${SLURM_JOB_ID:-}" ] && command -v scancel >/dev/null 2>&1; then
        log_warn "Issuing scancel for master job ${SLURM_JOB_ID}"
        scancel "${SLURM_JOB_ID}" >/dev/null 2>&1 || log_warn "scancel ${SLURM_JOB_ID} failed; exiting."
    fi
    exit 1
}

#
# Self-submit support
#
submit_self() {
    if [ -z "${DATASET:-}" ]; then
        echo "❌ DATASET is required before submitting to Slurm." >&2
        exit 1
    fi

    ensure_master_logging "${DATASET}" "prepare"

    # Build sbatch flags from pipeline config (account/partition/logs etc.)
    local script_path="${PIPELINE_ROOT}/bin/gatk_pipeline.sh"
    local stdout_template="${MASTER_LOG_DIR}/GATK_master_script_%j.output"
    local stderr_template="${MASTER_LOG_DIR}/GATK_master_script_%j.error"
    local sbatch_cmd=( sbatch
        -J GATK_master_script
        -A "${SLURM_ACCOUNT}"
        -p "${SLURM_PARTITION}"
        -N 1
        -n 1
        -c 1
        --mem=4G
        -t 336:00:00
        -o "${stdout_template}"
        -e "${stderr_template}"
        --export=ALL,PIPELINE_ROOT="${PIPELINE_ROOT}"
    )

    # Reconstruct arguments for inner execution
    local inner_args=( --as-sbatch )
    if [ -n "${DATASET:-}" ]; then
        inner_args+=( -d "${DATASET}" )
    fi
    if [ -n "${RDM_PATH:-}" ]; then
        inner_args+=( --rdm-path "${RDM_PATH}" )
    fi
    if [ -n "${REQUESTED_STEP:-}" ]; then
        inner_args+=( -s "${REQUESTED_STEP}" )
    fi
    if [ -n "${PIPELINE_INTERACTIVE_ACTION:-}" ]; then
        inner_args+=( --interactive-action "${PIPELINE_INTERACTIVE_ACTION}" )
    fi
    if [ "${PIPELINE_DRY_RUN:-false}" = "true" ]; then
        echo "[dry-run] Would execute: ${sbatch_cmd[*]} ${script_path} ${inner_args[*]}"
        echo "[dry-run] Master logs would be written to: ${MASTER_LOG_DIR}"
        return 0
    fi
    echo "Submitting master pipeline via sbatch..."
    echo "Master logs: ${MASTER_LOG_DIR}"
    "${sbatch_cmd[@]}" "${script_path}" "${inner_args[@]}"
}

_status_tail() {
    # Capture only the last non-empty line from a status function
    # to avoid mixing with logger output
    local out
    out="$("$@")" || true
    printf '%s\n' "${out}" | awk 'NF{p=$0}END{print p}'
}

wait_until_complete() {
    local label="$1"
    local checker_func="$2"
    local rdm_base="$3"
    local stall_count=0
    log_info "Waiting for ${label} to complete (poll ${poll_interval_seconds}s)..."
    while true; do
        local status
        status="$(_status_tail "${checker_func}" "${rdm_base}")"
        case "${status}" in
            "Complete")
                log_info "✅ ${label} complete."
                break
                ;;
            "Failed")
                log_error "❌ ${label} reported failure; aborting."
                return 1
                ;;
            "Partial"|"Incomplete"|"Not Started")
                stall_count=$((stall_count + 1))
                log_info "… ${label} status: ${status} @ $(date) (stall ${stall_count}/${master_stall_limit})"
                if [ "${stall_count}" -ge "${master_stall_limit}" ]; then
                    cancel_master_job "${label} not detected as running after ${stall_count} polls."
                fi
                sleep "${poll_interval_seconds}"
                ;;
            *)
                stall_count=$((stall_count + 1))
                log_warn "… ${label} status (unknown='${status}'), waiting… @ $(date) (stall ${stall_count}/${master_stall_limit})"
                if [ "${stall_count}" -ge "${master_stall_limit}" ]; then
                    cancel_master_job "${label} status unknown for ${stall_count} polls."
                fi
                sleep "${poll_interval_seconds}"
                ;;
        esac
    done
}

_step_state_dir() {
    local dataset_name="$1"
    pipeline_state_dir "${dataset_name}"
}

step1a_job_file() {
    local dataset_name="$1"
    echo "$(_step_state_dir "${dataset_name}")/step1a_job_id.txt"
}

_job_ids_active() {
    local job_file="$1"
    [ -f "${job_file}" ] || return 1
    while IFS= read -r job_id; do
        [ -n "${job_id}" ] || continue
        if squeue -h -j "${job_id}" >/dev/null 2>&1; then
            return 0
        fi
    done < "${job_file}"
    return 1
}

wait_for_step1a_outputs() {
    local dataset_name="$1"
    local rdm_base="$2"

    local state_dir
    state_dir="$(_step_state_dir "${dataset_name}")"
    mkdir -p "${state_dir}"
    local job_file="${state_dir}/step1a_job_id.txt"
    local heartbeat=0

    while true; do
        if _job_ids_active "${job_file}"; then
            if [ $((heartbeat % 5)) -eq 0 ]; then
                log_info "… Step 1A still running (job IDs from ${job_file})"
            fi
            heartbeat=$((heartbeat + 1))
            sleep "${poll_interval_seconds}"
            continue
        fi

        local s1a_status
        s1a_status="$(_status_tail check_step1a_status "${rdm_base}")"
        if [ "${s1a_status}" = "Complete" ]; then
            rm -f "${state_dir}/step1a_running.flag" "${state_dir}/step1a_failed.flag" 2>/dev/null || true
            log_info "Step 1A outputs verified as complete."
            return 0
        fi

        local missing
        missing="$(get_incomplete_samples "${rdm_base}")"
        local msg="Step 1A incomplete after jobs left queue."
        if [ -n "${missing}" ]; then
            msg+=" Missing samples: ${missing//$'\n'/, }"
        fi
        log_error "${msg}"
        printf '%s\n' "${msg}" > "${state_dir}/step1a_failed.flag"
        return 1
    done
}

run_master_loop() {
    local dataset_name="$1"
    local rdm_base="$2"
    local interactive_action="${PIPELINE_INTERACTIVE_ACTION:-}"

    ensure_master_logging "${dataset_name}"

    # Auto-approve to avoid interactive confirmations
    PIPELINE_AUTO_APPROVE="true"
    export PIPELINE_AUTO_APPROVE

    log_info "🧬 GATK master loop for dataset: ${dataset_name}"
    log_info "RDM base: ${rdm_base}"
    log_info "SLURM_JOB_ID: ${SLURM_JOB_ID:-none}"
    if [ -n "${interactive_action}" ]; then
        log_info "Interactive request detected: ${interactive_action}"
    fi

    case "${interactive_action}" in
        step1a)
        log_info "🚀 Running Step 1A per user request..."
            run_step1a "${dataset_name}" "${rdm_base}"
            wait_for_step1a_outputs "${dataset_name}" "${rdm_base}"
        unset PIPELINE_INTERACTIVE_ACTION
            ;;
        step1b)
            log_info "🚀 Running Step 1B per user request..."
            run_step1b "${dataset_name}" "${rdm_base}"
            wait_until_complete "Step 1B" check_step1b_status "${rdm_base}"
        unset PIPELINE_INTERACTIVE_ACTION
            ;;
        step1c)
            log_info "🚀 Running Step 1C per user request..."
            run_step1c "${dataset_name}" "${rdm_base}"
            wait_until_complete "Step 1C" check_step1c_status "${rdm_base}"
        unset PIPELINE_INTERACTIVE_ACTION
            ;;
        step1d)
            log_info "🚀 Running Step 1D per user request..."
            run_step1d "${dataset_name}" "${rdm_base}"
            wait_until_complete "Step 1D" check_step1d_status "${rdm_base}"
        unset PIPELINE_INTERACTIVE_ACTION
            ;;
        full)
            log_info "🚀 Running full pipeline per user request..."
            run_complete_pipeline "${dataset_name}" "${rdm_base}"
        unset PIPELINE_INTERACTIVE_ACTION
            return
            ;;
    esac

    # Step 1A
    local s1a_status
    s1a_status="$(_status_tail check_step1a_status "${rdm_base}")"
    if [ "${s1a_status}" != "Complete" ]; then
        local job_file
        job_file="$(step1a_job_file "${dataset_name}")"
        if [ -f "${job_file}" ] && _job_ids_active "${job_file}"; then
            log_info "Step 1A already submitted (job IDs in ${job_file}); waiting for completion."
        else
            log_info "🚀 Submitting Step 1A…"
            if ! run_step1a "${dataset_name}" "${rdm_base}"; then
                log_error "❌ Step 1A submission failed; aborting master loop."
                exit 1
            fi
        fi
        if ! wait_for_step1a_outputs "${dataset_name}" "${rdm_base}"; then
            exit 1
        fi
    else
        log_info "✅ Step 1A already complete."
    fi

    # Step 1B
    local s1b_status
    s1b_status="$(_status_tail check_step1b_status "${rdm_base}")"
    if [ "${s1b_status}" != "Complete" ]; then
        log_info "🚀 Submitting Step 1B…"
        if ! run_step1b "${dataset_name}" "${rdm_base}"; then
            log_error "❌ Step 1B submission failed; aborting master loop."
            exit 1
        fi
        if ! wait_until_complete "Step 1B" check_step1b_status "${rdm_base}"; then
            exit 1
        fi
    else
        log_info "✅ Step 1B already complete."
    fi

    # Step 1C
    local s1c_status
    s1c_status="$(_status_tail check_step1c_status "${rdm_base}")"
    if [ "${s1c_status}" != "Complete" ]; then
        log_info "🚀 Submitting Step 1C…"
        if ! run_step1c "${dataset_name}" "${rdm_base}"; then
            log_error "❌ Step 1C submission failed; aborting master loop."
            exit 1
        fi
        if ! wait_until_complete "Step 1C" check_step1c_status "${rdm_base}"; then
            exit 1
        fi
    else
        log_info "✅ Step 1C already complete."
    fi

    # Step 1D
    local s1d_status
    s1d_status="$(_status_tail check_step1d_status "${rdm_base}")"
    if [ "${s1d_status}" != "Complete" ]; then
        log_info "🚀 Submitting Step 1D…"
        if ! run_step1d "${dataset_name}" "${rdm_base}"; then
            log_error "❌ Step 1D submission failed; aborting master loop."
            exit 1
        fi
        if ! wait_until_complete "Step 1D" check_step1d_status "${rdm_base}"; then
            exit 1
        fi
    else
        log_info "✅ Step 1D already complete."
    fi

    log_info "🎉 All steps complete for ${dataset_name}."
}

print_usage() {
    cat <<'USAGE'
Usage: bin/gatk_pipeline.sh [options]

Options:
  -d, --dataset NAME        Dataset name (required for non-interactive mode)
      --rdm-path PATH       Override RDM base path
  -s, --step STEP           Step to run (1a,1b,1c,1d,full,auto) [default: auto]
  -i, --interactive         Force interactive mode
      --dry-run             Show actions without submitting jobs
  -q, --quiet               Suppress banners/prompts where possible
      --submit              Submit this script to Slurm (self-submit mode)
      --as-sbatch           Internal flag: indicates running inside Slurm allocation
  -h, --help                Show this message and exit
USAGE
}

run_cli_mode() {
    local dataset_name="$DATASET"
    local rdm_base="${RDM_PATH:-${RDM_DATASETS_PATH}/${DATASET}}"

    ensure_master_logging "${dataset_name}"

    if [ ! -d "${rdm_base}" ]; then
        log_error "❌ RDM path not found: ${rdm_base}"
        exit 1
    fi

    export PIPELINE_DRY_RUN PIPELINE_QUIET

    case "$REQUESTED_STEP" in
        1a)
            run_step1a "${dataset_name}" "${rdm_base}" ;;
        1b)
            run_step1b "${dataset_name}" "${rdm_base}" ;;
        1c)
            run_step1c "${dataset_name}" "${rdm_base}" ;;
        1d)
            run_step1d "${dataset_name}" "${rdm_base}" ;;
        full)
            run_complete_pipeline "${dataset_name}" "${rdm_base}" ;;
        auto)
            PIPELINE_AUTO_APPROVE="true"
            export PIPELINE_AUTO_APPROVE
            local dataset_info="${dataset_name}|${rdm_base}"
            local status=$(analyze_pipeline_status "${dataset_info}")
            print_user_report "${dataset_name}" "${rdm_base}" "${status}"
            route_pipeline_execution "${dataset_info}" "${status}"
            ;;
        *)
            log_error "❌ Unknown step: ${REQUESTED_STEP}"
            exit 1;;
    esac
}

run_interactive_mode() {
    $PIPELINE_QUIET || show_welcome
    local dataset_info
    dataset_info=$(get_dataset_info)
    local status
    status=$(analyze_pipeline_status "$dataset_info")

    # Extract dataset name for submission
    DATASET="$(printf '%s' "$dataset_info" | cut -d'|' -f1)"
    local RDM_BASE_INTERACTIVE
    RDM_BASE_INTERACTIVE="$(printf '%s' "$dataset_info" | cut -d'|' -f2)"
    REQUESTED_STEP="${REQUESTED_STEP:-auto}"

    # Resolve actions interactively before submitting
    if [ -z "${SLURM_JOB_ID:-}" ]; then
        ensure_master_logging "${DATASET}" "prepare"
        print_user_report "${DATASET}" "${RDM_BASE_INTERACTIVE}" "${status}"
        route_pipeline_execution "$dataset_info" "$status" "interactive"
        if [ -n "${PIPELINE_INTERACTIVE_ACTION:-}" ]; then
            echo "➡️  Launching selected action (${PIPELINE_INTERACTIVE_ACTION}) via Slurm..."
            submit_self
            exit 0
        else
            echo "✅ No pending actions selected; exiting interactive session."
            exit 0
        fi
    fi

    ensure_master_logging "${DATASET}"
    route_pipeline_execution "$dataset_info" "$status"
}

DATASET=""
RDM_PATH=""
REQUESTED_STEP="auto"
FORCE_INTERACTIVE=false
PIPELINE_DRY_RUN="false"
PIPELINE_QUIET="false"
SUBMIT_MODE=false
AS_SBATCH=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        -d|--dataset)
            DATASET="$2"; shift 2;;
        --rdm-path)
            RDM_PATH="$2"; shift 2;;
        -s|--step)
            REQUESTED_STEP=$(echo "$2" | tr '[:upper:]' '[:lower:]'); shift 2;;
        -i|--interactive)
            FORCE_INTERACTIVE=true; shift;;
        --dry-run)
            PIPELINE_DRY_RUN="true"; shift;;
        -q|--quiet)
            PIPELINE_QUIET="true"; shift;;
        --submit)
            SUBMIT_MODE=true; shift;;
        --as-sbatch)
            AS_SBATCH=true; shift;;
        --interactive-action)
            PIPELINE_INTERACTIVE_ACTION="$2"; shift 2;;
        -h|--help)
            print_usage; exit 0;;
        *)
            echo "Unknown option: $1" >&2
            print_usage
            exit 1;;
    esac
done

export PIPELINE_DRY_RUN PIPELINE_QUIET

# If explicitly asked to submit, do so and exit
if $SUBMIT_MODE; then
    submit_self
    exit 0
fi

# If running inside Slurm (via --as-sbatch or allocated job), execute non-interactive paths
if $AS_SBATCH || [ -n "${SLURM_JOB_ID:-}" ]; then
    if ! $FORCE_INTERACTIVE && [ -n "$DATASET" ] && [ "$REQUESTED_STEP" = "auto" ]; then
        # Resolve RDM base path (same logic as CLI mode)
        rdm_base="${RDM_PATH:-${RDM_DATASETS_PATH}/${DATASET}}"
        ensure_master_logging "${DATASET}"
        if [ ! -d "${rdm_base}" ]; then
            log_error "❌ RDM path not found: ${rdm_base}"
            exit 1
        fi
        run_master_loop "${DATASET}" "${rdm_base}"
    elif $FORCE_INTERACTIVE || [ -z "$DATASET" ]; then
        run_interactive_mode
    else
        run_cli_mode
    fi
    exit 0
fi

if ! $FORCE_INTERACTIVE && [ -n "$DATASET" ] && [ "$REQUESTED_STEP" = "auto" ]; then
    # Resolve RDM base path (same logic as CLI mode)
    rdm_base="${RDM_PATH:-${RDM_DATASETS_PATH}/${DATASET}}"
    ensure_master_logging "${DATASET}"
    if [ ! -d "${rdm_base}" ]; then
        log_error "❌ RDM path not found: ${rdm_base}"
        exit 1
    fi
    run_master_loop "${DATASET}" "${rdm_base}"
elif $FORCE_INTERACTIVE || [ -z "$DATASET" ]; then
    run_interactive_mode
else
    run_cli_mode
fi
