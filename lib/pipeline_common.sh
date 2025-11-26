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

# Function to run Step 1A (calls modular backend)
run_step1a() {
    local dataset_name="$1"
    local rdm_base_path="$2"
    local prev_log_override="${PIPELINE_LOG_DIR_OVERRIDE:-}"
    
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
