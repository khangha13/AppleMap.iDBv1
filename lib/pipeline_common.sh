#!/bin/bash
# =============================================================================
# COMMON PIPELINE FUNCTIONS FOR GATK PIPELINE
# =============================================================================
# Common functions used across all pipeline steps

# Function to get dataset information from user
get_dataset_info() {
    {
        echo "📋 DATASET CONFIGURATION"
        echo "========================="
        echo ""
    } >&2
    
    # Get dataset name
    local dataset_name
    while true; do
        read -p "Please state your dataset: " dataset_name
        if [ -n "$dataset_name" ]; then
            break
        else
            echo "❌ Dataset name cannot be empty. Please try again." >&2
        fi
    done
    
    # Infer path
    local rdm_base_path="${RDM_DATASETS_PATH}/${dataset_name}"
    
    {
        echo ""
        echo "🔍 ANALYZING DATASET..."
        echo "Dataset: $dataset_name"
        echo "Path: $rdm_base_path"
    } >&2
    
    # Validate path exists
    if [ ! -d "$rdm_base_path" ]; then
        echo "❌ Error: Path does not exist: $rdm_base_path" >&2
        echo "Please ensure the directory exists and contains the required structure." >&2
        exit 1
    fi
    
    # Return dataset info
    echo "$dataset_name|$rdm_base_path"
}

# Function to analyze pipeline status
analyze_pipeline_status() {
    local dataset_info="$1"
    local dataset_name=$(echo "$dataset_info" | cut -d'|' -f1)
    local rdm_base_path=$(echo "$dataset_info" | cut -d'|' -f2)
    
    {
        echo ""
        echo "🔍 PIPELINE STATUS ANALYSIS"
        echo "============================"
    } >&2
    
    local step1a_status=$(check_step1a_status "$rdm_base_path")
    echo "Step 1A (Per-sample): $step1a_status" >&2
    local step1b_status=$(check_step1b_status "$rdm_base_path")
    echo "Step 1B (CombineGVCFs): $step1b_status" >&2
    local step1c_status=$(check_step1c_status "$rdm_base_path")
    echo "Step 1C (Beagle Imputation): $step1c_status" >&2
    local step1d_status=$(check_step1d_status "$rdm_base_path")
    echo "Step 1D (QC Analysis): $step1d_status" >&2

    if [[ "$step1a_status" != "Complete" ]]; then
        echo "Overall Status: 🔄 Step 1A Needed" >&2
        echo "step1a_needed"
    elif [[ "$step1b_status" != "Complete" ]]; then
        echo "Overall Status: ⚠️  Step 1B Needed" >&2
        echo "step1b_needed"
    elif [[ "$step1c_status" != "Complete" ]]; then
        echo "Overall Status: ⚠️  Step 1C Needed" >&2
        echo "step1c_needed"
    elif [[ "$step1d_status" != "Complete" ]]; then
        echo "Overall Status: ⚠️  Step 1D Needed" >&2
        echo "step1d_needed"
    else
        echo "Overall Status: ✅ Complete" >&2
        echo "complete"
    fi
}

# Function to handle complete pipeline
handle_complete_pipeline() {
    local dataset_name="$1"
    local rdm_base_path="$2"
    local mode="${3:-auto}"
    
    echo "🎉 Congratulations! Your pipeline is complete."
    echo ""
    echo "Results available at:"
    echo "• Individual VCFs: $rdm_base_path/5.Individual_VCF/"
    echo "• Consolidated VCFs: $rdm_base_path/7.Consolidated_VCF/"
    echo ""
    
    if [ "$mode" = "interactive" ]; then
        echo "✅ No actions required. Exiting interactive session."
    elif confirm "Would you like to rerun the entire pipeline?"; then
        echo "🔄 Proceeding with full pipeline rerun..."
        run_complete_pipeline "$dataset_name" "$rdm_base_path"
    else
        echo "✅ Pipeline execution complete. Goodbye!"
        exit 0
    fi
}

# Function to handle Step 1B needed
handle_step1b_needed() {
    local dataset_name="$1"
    local rdm_base_path="$2"
    local mode="${3:-auto}"
    
    echo "Step 1A is complete, but Step 1B (CombineGVCFs) needs to be run."
    echo ""
    
    if [ "$mode" = "interactive" ]; then
        if confirm "Would you like to run Step 1B now?"; then
            PIPELINE_INTERACTIVE_ACTION="step1b"
            PIPELINE_INTERACTIVE_DATASET="$dataset_name"
            PIPELINE_INTERACTIVE_RDM="$rdm_base_path"
            export PIPELINE_INTERACTIVE_ACTION PIPELINE_INTERACTIVE_DATASET PIPELINE_INTERACTIVE_RDM
            echo "🚀 Step 1B queued for submission."
        else
            echo "✅ Step 1B execution skipped. You can run it later."
        fi
        return
    fi
    
    if confirm "Would you like to run Step 1B now?"; then
        echo "🚀 Running Step 1B..."
        run_step1b "$dataset_name" "$rdm_base_path"
    else
        echo "✅ Step 1B execution cancelled. You can run it later."
        exit 0
    fi
}

# Function to handle Step 1A needed
handle_step1a_needed() {
    local dataset_name="$1"
    local rdm_base_path="$2"
    local mode="${3:-auto}"
    
    echo "Step 1A (Per-sample variant calling) needs to be run."
    echo ""
    
    if [ "$mode" = "interactive" ]; then
        if confirm "Would you like to run Step 1A now?"; then
            PIPELINE_INTERACTIVE_ACTION="step1a"
            PIPELINE_INTERACTIVE_DATASET="$dataset_name"
            PIPELINE_INTERACTIVE_RDM="$rdm_base_path"
            export PIPELINE_INTERACTIVE_ACTION PIPELINE_INTERACTIVE_DATASET PIPELINE_INTERACTIVE_RDM
            echo "🚀 Step 1A queued for submission."
        else
            echo "✅ Step 1A execution skipped. You can run it later."
        fi
        return
    fi
    
    if confirm "Would you like to run Step 1A now?"; then
        echo "🚀 Running Step 1A..."
        run_step1a "$dataset_name" "$rdm_base_path"
    else
        echo "✅ Step 1A execution cancelled. You can run it later."
        exit 0
    fi
}

# Function to handle unknown status
handle_step1c_needed() {
    local dataset_name="$1"
    local rdm_base_path="$2"
    local mode="${3:-auto}"

    echo "Steps 1A and 1B are complete, but Step 1C (Beagle imputation) needs to be run."
    echo ""

    if [ "$mode" = "interactive" ]; then
        if confirm "Would you like to run Step 1C now?"; then
            PIPELINE_INTERACTIVE_ACTION="step1c"
            PIPELINE_INTERACTIVE_DATASET="$dataset_name"
            PIPELINE_INTERACTIVE_RDM="$rdm_base_path"
            export PIPELINE_INTERACTIVE_ACTION PIPELINE_INTERACTIVE_DATASET PIPELINE_INTERACTIVE_RDM
            echo "🚀 Step 1C queued for submission."
        else
            echo "✅ Step 1C execution skipped. You can run it later."
        fi
        return
    fi

    if confirm "Would you like to run Step 1C now?"; then
        echo "🚀 Running Step 1C..."
        run_step1c "$dataset_name" "$rdm_base_path"
    else
        echo "✅ Step 1C execution cancelled. You can run it later."
        exit 0
    fi
}

handle_step1d_needed() {
    local dataset_name="$1"
    local rdm_base_path="$2"
    local mode="${3:-auto}"

    echo "Quality-control analysis (Step 1D) has not been run yet."
    echo ""

    if [ "$mode" = "interactive" ]; then
        if confirm "Would you like to run Step 1D now?"; then
            PIPELINE_INTERACTIVE_ACTION="step1d"
            PIPELINE_INTERACTIVE_DATASET="$dataset_name"
            PIPELINE_INTERACTIVE_RDM="$rdm_base_path"
            export PIPELINE_INTERACTIVE_ACTION PIPELINE_INTERACTIVE_DATASET PIPELINE_INTERACTIVE_RDM
            echo "🚀 Step 1D queued for submission."
        else
            echo "✅ Step 1D execution skipped. You can run it later."
        fi
        return
    fi

    if confirm "Would you like to run Step 1D now?"; then
        echo "🚀 Running Step 1D..."
        run_step1d "$dataset_name" "$rdm_base_path"
    else
        echo "✅ Step 1D execution cancelled. You can run it later."
        exit 0
    fi
}

handle_unknown_status() {
    local dataset_name="$1"
    local rdm_base_path="$2"
    local mode="${3:-auto}"
    
    echo "❓ Unable to determine pipeline status."
    echo "This might be a new dataset or the directory structure is incomplete."
    echo ""
    
    if [ "$mode" = "interactive" ]; then
        if confirm "Would you like to run the complete pipeline from the beginning?"; then
            PIPELINE_INTERACTIVE_ACTION="full"
            PIPELINE_INTERACTIVE_DATASET="$dataset_name"
            PIPELINE_INTERACTIVE_RDM="$rdm_base_path"
            export PIPELINE_INTERACTIVE_ACTION PIPELINE_INTERACTIVE_DATASET PIPELINE_INTERACTIVE_RDM
            echo "🚀 Full pipeline queued for submission."
        else
            echo "✅ Pipeline execution cancelled."
        fi
        return
    fi
    
    if confirm "Would you like to run the complete pipeline from the beginning?"; then
        echo "🚀 Running complete pipeline..."
        run_complete_pipeline "$dataset_name" "$rdm_base_path"
    else
        echo "✅ Pipeline execution cancelled."
        exit 0
    fi
}

# Simple confirmation function
confirm() {
    local message="$1"
    local response
    if [ "${PIPELINE_AUTO_APPROVE:-false}" = "true" ]; then
        echo "${message} [auto-yes]"
        return 0
    fi
    
    while true; do
        read -p "$message [y/n]: " response
        case $response in
            [Yy]* ) return 0;;
            [Nn]* ) return 1;;
            * ) echo "Please answer yes or no.";;
        esac
    done
}

# -----------------------------------------------------------------------------
# Shared reference preparation (moved to master preflight)
# -----------------------------------------------------------------------------
pipeline_state_dir() {
    local dataset_name="$1"
    if [ -n "${MASTER_LOG_DIR:-}" ]; then
        echo "${MASTER_LOG_DIR}"
    else
        echo "${LOG_BASE_PATH%/}/${dataset_name}"
    fi
}

shared_reference_base_path() {
    local dataset_name="$1"
    echo "${SCRATCH_BASE_PATH%/}/${dataset_name}_shared"
}

_require_source_file() {
    local path="$1"
    local label="$2"
    if [ ! -f "${path}" ]; then
        echo "${label} (${path})"
    fi
}

validate_shared_reference_manifest() {
    local dataset_name="$1"
    local reference_genome="${2:-$(get_reference_fasta)}"
    local known_sites="${3:-$(get_known_sites_vcf)}"
    local adapter_file="${4:-$(get_adapter_fasta)}"

    local shared_base
    shared_base="$(shared_reference_base_path "${dataset_name}")"
    local shared_ref_dir="${shared_base}/Reference_genome"
    local shared_known_dir="${shared_base}/Known_sites"
    local shared_adapter_dir="${shared_base}/Adapter_file"

    local ref_basename known_sites_basename adapter_basename
    ref_basename="$(basename "${reference_genome}")"
    known_sites_basename="$(basename "${known_sites}")"
    adapter_basename="$(basename "${adapter_file}")"
    local ref_dict_basename="${reference_genome%.*}.dict"
    ref_dict_basename="$(basename "${ref_dict_basename}")"

    local missing=()
    [ -f "${shared_ref_dir}/${ref_basename}" ] || missing+=("reference FASTA (${shared_ref_dir}/${ref_basename})")
    [ -f "${shared_ref_dir}/${ref_basename}.fai" ] || missing+=("reference FASTA index (.fai)")
    [ -f "${shared_ref_dir}/${ref_dict_basename}" ] || missing+=("reference dictionary (.dict)")
    for ext in amb ann bwt pac sa; do
        [ -f "${shared_ref_dir}/${ref_basename}.${ext}" ] || missing+=("BWA index .${ext}")
    done
    [ -f "${shared_known_dir}/${known_sites_basename}" ] || missing+=("known sites VCF (${shared_known_dir}/${known_sites_basename})")
    local has_known_index="false"
    for idx_ext in tbi idx; do
        if [ -f "${shared_known_dir}/${known_sites_basename}.${idx_ext}" ]; then
            has_known_index="true"
            break
        fi
    done
    if [ "${has_known_index}" != "true" ]; then
        missing+=("known sites index (.tbi or .idx)")
    fi
    [ -f "${shared_adapter_dir}/${adapter_basename}" ] || missing+=("adapter FASTA (${shared_adapter_dir}/${adapter_basename})")

    if [ ${#missing[@]} -gt 0 ]; then
        log_error "Shared reference assets are incomplete for dataset ${dataset_name}: ${missing[*]}"
        return 1
    fi

    return 0
}

ensure_shared_references_ready() {
    local dataset_name="$1"
    local reference_genome="${2:-$(get_reference_fasta)}"
    local known_sites="${3:-$(get_known_sites_vcf)}"
    local adapter_file="${4:-$(get_adapter_fasta)}"

    if ! validate_shared_reference_manifest "${dataset_name}" "${reference_genome}" "${known_sites}" "${adapter_file}"; then
        error_exit "Shared reference assets missing for dataset ${dataset_name}. Please rerun the master script to stage references."
    fi
}

prepare_shared_reference_assets() {
    local dataset_name="$1"
    local reference_genome="${2:-$(get_reference_fasta)}"
    local known_sites="${3:-$(get_known_sites_vcf)}"
    local adapter_file="${4:-$(get_adapter_fasta)}"

    if [ -z "${dataset_name:-}" ]; then
        error_exit "prepare_shared_reference_assets requires a dataset name."
    fi

    local shared_base
    shared_base="$(shared_reference_base_path "${dataset_name}")"
    local shared_ref_dir="${shared_base}/Reference_genome"
    local shared_known_dir="${shared_base}/Known_sites"
    local shared_adapter_dir="${shared_base}/Adapter_file"
    local marker_file="${shared_base}/.initialized"

    if validate_shared_reference_manifest "${dataset_name}" "${reference_genome}" "${known_sites}" "${adapter_file}"; then
        log_info "Shared reference assets already prepared at ${shared_base}"
        return 0
    fi

    if [ "${PIPELINE_DRY_RUN:-false}" = "true" ]; then
        log_info "[dry-run] Would stage shared reference assets to ${shared_base}"
        return 0
    fi

    local missing_sources=()
    local source_missing
    source_missing="$(_require_source_file "${reference_genome}" "reference FASTA")"; [ -n "${source_missing}" ] && missing_sources+=("${source_missing}")
    source_missing="$(_require_source_file "${reference_genome}.fai" "reference FASTA index (.fai)")"; [ -n "${source_missing}" ] && missing_sources+=("${source_missing}")
    local ref_dict="${reference_genome%.*}.dict"
    source_missing="$(_require_source_file "${ref_dict}" "reference dictionary (.dict)")"; [ -n "${source_missing}" ] && missing_sources+=("${source_missing}")
    for ext in amb ann bwt pac sa; do
        source_missing="$(_require_source_file "${reference_genome}.${ext}" "BWA index .${ext}")"; [ -n "${source_missing}" ] && missing_sources+=("${source_missing}")
    done
    source_missing="$(_require_source_file "${known_sites}" "known sites VCF")"; [ -n "${source_missing}" ] && missing_sources+=("${source_missing}")
    local known_sites_index=""
    if [ -f "${known_sites}.tbi" ]; then
        known_sites_index="${known_sites}.tbi"
    elif [ -f "${known_sites}.idx" ]; then
        known_sites_index="${known_sites}.idx"
    fi
    if [ -z "${known_sites_index}" ]; then
        missing_sources+=("known sites index (.tbi or .idx)")
    fi
    source_missing="$(_require_source_file "${adapter_file}" "adapter FASTA")"; [ -n "${source_missing}" ] && missing_sources+=("${source_missing}")

    if [ ${#missing_sources[@]} -gt 0 ]; then
        error_exit "Cannot stage shared references; missing source files: ${missing_sources[*]}"
    fi

    log_info "Staging shared reference assets to ${shared_base}"
    mkdir -p "${shared_ref_dir}" "${shared_known_dir}" "${shared_adapter_dir}"
    rm -f "${marker_file}"

    rsync -rhPt "${reference_genome}" "${shared_ref_dir}/" || error_exit "Failed to copy reference FASTA to shared location"
    rsync -rhPt "${reference_genome}.fai" "${shared_ref_dir}/" || error_exit "Failed to copy reference FASTA index (.fai) to shared location"
    rsync -rhPt "${ref_dict}" "${shared_ref_dir}/" || error_exit "Failed to copy reference dictionary (.dict) to shared location"
    for ext in amb ann bwt pac sa alt; do
        if [ -f "${reference_genome}.${ext}" ]; then
            rsync -rhPt "${reference_genome}.${ext}" "${shared_ref_dir}/" || error_exit "Failed to copy BWA index .${ext} to shared location"
        fi
    done

    rsync -rhPt "${known_sites}" "${shared_known_dir}/" || error_exit "Failed to copy known sites VCF to shared location"
    if [ -n "${known_sites_index}" ]; then
        rsync -rhPt "${known_sites_index}" "${shared_known_dir}/" || error_exit "Failed to copy known sites index to shared location"
    fi

    rsync -rhPt "${adapter_file}" "${shared_adapter_dir}/" || error_exit "Failed to copy adapter FASTA to shared location"

    if ! validate_shared_reference_manifest "${dataset_name}" "${reference_genome}" "${known_sites}" "${adapter_file}"; then
        error_exit "Shared reference assets still incomplete after staging to ${shared_base}"
    fi

    touch "${marker_file}"
    log_info "Shared reference assets prepared at ${shared_base}"
}

# Function to run Step 1A (calls modular backend)
run_step1a() {
    local dataset_name="$1"
    local rdm_base_path="$2"
    local prev_log_override="${PIPELINE_LOG_DIR_OVERRIDE:-}"
    
    prepare_shared_reference_assets "${dataset_name}"

    if [ "${PIPELINE_DRY_RUN:-false}" = "true" ]; then
        echo "[dry-run] Step 1A would run for ${dataset_name} (${rdm_base_path})."
        return 0
    fi
    if [ -z "${PIPELINE_LOG_DIR_OVERRIDE:-}" ] && [ -n "${MASTER_LOG_DIR:-}" ]; then
        export PIPELINE_LOG_DIR_OVERRIDE="${MASTER_LOG_DIR}"
    fi
    bash "${PIPELINE_ROOT}/wrappers/sbatch/step1a_submit.sh" "$dataset_name" "$rdm_base_path"
    local rc=$?
    if [ -n "${prev_log_override}" ]; then
        PIPELINE_LOG_DIR_OVERRIDE="${prev_log_override}"
    else
        unset PIPELINE_LOG_DIR_OVERRIDE 2>/dev/null || true
    fi
    if [ $rc -ne 0 ]; then
        log_error "Step 1A submission wrapper exited with code ${rc}."
    fi
    return $rc
}

# Function to run Step 1B (calls modular backend)
run_step1b() {
    local dataset_name="$1"
    local rdm_base_path="$2"
    local prev_log_override="${PIPELINE_LOG_DIR_OVERRIDE:-}"

    prepare_shared_reference_assets "${dataset_name}"
    
    if [ "${PIPELINE_DRY_RUN:-false}" = "true" ]; then
        echo "[dry-run] Step 1B would run for ${dataset_name}."
        return 0
    fi
    if [ -z "${PIPELINE_LOG_DIR_OVERRIDE:-}" ] && [ -n "${MASTER_LOG_DIR:-}" ]; then
        export PIPELINE_LOG_DIR_OVERRIDE="${MASTER_LOG_DIR}"
    fi
    bash "${PIPELINE_ROOT}/wrappers/sbatch/step1b_submit.sh" "$dataset_name" "$rdm_base_path"
    local rc=$?
    if [ -n "${prev_log_override}" ]; then
        PIPELINE_LOG_DIR_OVERRIDE="${prev_log_override}"
    else
        unset PIPELINE_LOG_DIR_OVERRIDE 2>/dev/null || true
    fi
    if [ $rc -ne 0 ]; then
        log_error "Step 1B submission wrapper exited with code ${rc}."
    fi
    return $rc
}

# Function to run complete pipeline
run_step1c() {
    local dataset_name="$1"
    local rdm_base_path="$2"
    local prev_log_override="${PIPELINE_LOG_DIR_OVERRIDE:-}"
    if [ "${PIPELINE_DRY_RUN:-false}" = "true" ]; then
        echo "[dry-run] Step 1C would run for ${dataset_name}."
        return 0
    fi
    if [ -z "${PIPELINE_LOG_DIR_OVERRIDE:-}" ] && [ -n "${MASTER_LOG_DIR:-}" ]; then
        export PIPELINE_LOG_DIR_OVERRIDE="${MASTER_LOG_DIR}"
    fi
    bash "${PIPELINE_ROOT}/wrappers/sbatch/step1c_submit.sh" "$dataset_name" "$rdm_base_path"
    local rc=$?
    if [ -n "${prev_log_override}" ]; then
        PIPELINE_LOG_DIR_OVERRIDE="${prev_log_override}"
    else
        unset PIPELINE_LOG_DIR_OVERRIDE 2>/dev/null || true
    fi
    if [ $rc -ne 0 ]; then
        log_error "Step 1C submission wrapper exited with code ${rc}."
    fi
    return $rc
}

run_step1d() {
    local dataset_name="$1"
    local rdm_base_path="$2"
    local prev_log_override="${PIPELINE_LOG_DIR_OVERRIDE:-}"
    if [ "${PIPELINE_DRY_RUN:-false}" = "true" ]; then
        echo "[dry-run] Step 1D would run for ${dataset_name}."
        return 0
    fi
    if [ -z "${PIPELINE_LOG_DIR_OVERRIDE:-}" ] && [ -n "${MASTER_LOG_DIR:-}" ]; then
        export PIPELINE_LOG_DIR_OVERRIDE="${MASTER_LOG_DIR}"
    fi
    bash "${PIPELINE_ROOT}/wrappers/sbatch/step1d_submit.sh" "$dataset_name" "$rdm_base_path"
    local rc=$?
    if [ -n "${prev_log_override}" ]; then
        PIPELINE_LOG_DIR_OVERRIDE="${prev_log_override}"
    else
        unset PIPELINE_LOG_DIR_OVERRIDE 2>/dev/null || true
    fi
    if [ $rc -ne 0 ]; then
        log_error "Step 1D submission wrapper exited with code ${rc}."
    fi
    return $rc
}

run_complete_pipeline() {
    local dataset_name="$1"
    local rdm_base_path="$2"
    
    echo "🔄 Running complete pipeline (Step 1A → Step 1B → Step 1C → Step 1D)..."
    
    run_step1a "$dataset_name" "$rdm_base_path" || return 1
    run_step1b "$dataset_name" "$rdm_base_path" || return 1
    run_step1c "$dataset_name" "$rdm_base_path" || return 1
    run_step1d "$dataset_name" "$rdm_base_path" || return 1
}

# Function to show welcome message
show_welcome() {
    clear
    echo ""
    echo "╔══════════════════════════════════════════════════════════════════════════════════╗"
    echo "║                                                                                  ║"
    echo "║                    🧬 GATK VARIANT CALLING PIPELINE 🧬                        ║"
    echo "║                                                                                  ║"
    echo "║  Welcome! This pipeline will guide you through the complete GATK workflow:      ║"
    echo "║  • Step 1A: Per-sample variant calling (HaplotypeCaller)                       ║"
    echo "║  • Step 1B: Combine GVCFs (GenomicsDBImport + GenotypeGVCFs)                   ║"
    echo "║  • Step 1C: Beagle imputation                                                   ║"
    echo "║  • Step 1D: Post-imputation QC & analytics                                      ║"
    echo "║                                                                                  ║"
    echo "║  The pipeline will automatically detect what needs to be done and guide you.  ║"
    echo "║                                                                                  ║"
    echo "╚══════════════════════════════════════════════════════════════════════════════════╝"
    echo ""
}

# Function to print a short user summary
print_user_report() {
    local dataset_name="$1"
    local rdm_base_path="$2"
    local pipeline_status="$3"

    local status_message=""
    case "$pipeline_status" in
        step1a_needed) status_message="Step 1A pending";;
        step1b_needed) status_message="Step 1B pending";;
        step1c_needed) status_message="Step 1C pending";;
        step1d_needed) status_message="Step 1D pending";;
        complete)      status_message="All steps complete";;
        *)             status_message="Status unknown";;
    esac

    cat <<EOF
📊 DATASET SUMMARY
------------------
Name : ${dataset_name}
Path : ${rdm_base_path}
State: ${status_message}
EOF
}

# Function to route pipeline execution based on status
route_pipeline_execution() {
    local dataset_info="$1"
    local pipeline_status="$2"
    local mode="${3:-auto}"
    local dataset_name=$(echo "$dataset_info" | cut -d'|' -f1)
    local rdm_base_path=$(echo "$dataset_info" | cut -d'|' -f2)
    
    echo ""
    echo "🚀 PIPELINE EXECUTION"
    echo "====================="
    
    case "$pipeline_status" in
        "complete")
            handle_complete_pipeline "$dataset_name" "$rdm_base_path" "$mode"
            ;;
        "step1b_needed")
            handle_step1b_needed "$dataset_name" "$rdm_base_path" "$mode"
            ;;
        "step1a_needed")
            handle_step1a_needed "$dataset_name" "$rdm_base_path" "$mode"
            ;;
        "step1c_needed")
            handle_step1c_needed "$dataset_name" "$rdm_base_path" "$mode"
            ;;
        "step1d_needed")
            handle_step1d_needed "$dataset_name" "$rdm_base_path" "$mode"
            ;;
        "unknown")
            handle_unknown_status "$dataset_name" "$rdm_base_path" "$mode"
            ;;
    esac
}
