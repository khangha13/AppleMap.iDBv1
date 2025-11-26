#!/bin/bash

# =============================================================================
# BAM MOSDEPTH ANALYSIS PIPELINE
# =============================================================================
# 
# DESCRIPTION:
#   This script performs coverage analysis on BAM files using mosdepth,
#   then merges the results and generates visualization plots for all samples
#   in a dataset. It handles the complete workflow from BAM files to final plots.
#
# USAGE:
#   Normal run (submits SLURM job for all samples):
#     bash BAM_mosdepth_pipeline.sh
#
#   Dry run (preview without submitting):
#     bash BAM_mosdepth_pipeline.sh --dry-run
#
#   Test run with a single sample:
#     bash BAM_mosdepth_pipeline.sh --sample <SAMPLE_NAME>
#     bash BAM_mosdepth_pipeline.sh --sample <SAMPLE_NAME> --dry-run
#
# WHAT IT DOES:
#   1. Prompts for dataset name and validates the path
#   2. Scans for recalibrated BAM files (*_recal.bam) in the dataset's 4.BAM directory
#      - Normal mode: Finds ALL samples in the dataset
#      - Test mode (--sample): Finds only the specified sample
#   3. Creates a sample list file with recalibrated BAM file path(s)
#   4. Generates a SLURM script that:
#      - Copies all BAM files to TMPDIR (local SSD for fast processing)
#      - Runs mosdepth on each sample (2 threads per sample) with fast mode
#      - Runs in fast mode (-n --fast-mode) with 250bp bins (--by 250)
#      - Merges all coverage files from 250bp bins
#      - Generates visualization plots (faceted and individual chromosome plots)
#      - Copies results back to /scratch/user/uqpha1/mosdepth/
#   5. Submits the SLURM job or previews it (dry-run mode)
#   6. Adjusts resources automatically: lighter for test runs, full for all samples
#
# OUTPUT:
#   - Sample list file: bam_list_<DATASET>_<TIMESTAMP>.txt
#   - SLURM script: mosdepth_<DATASET>_<TIMESTAMP>.sh
#   - Job results in: /scratch/user/uqpha1/mosdepth/<DATASET>_merged/
#     * merged_coverage.bed.gz - merged coverage data
#     * merged_coverage_faceted.png - all chromosomes in one plot
#     * individual_chromosomes_combined.png - individual chromosome plots
#
# REQUIREMENTS:
#   - Access to /QRISdata/Q8367/WGS_Reference_Panel/
#   - Access to /scratch/user/uqpha1/
#   - SLURM cluster (e.g., Bunya at UQ)
#   - Dataset must have BAM files in: /QRISdata/Q8367/WGS_Reference_Panel/<DATASET>/4.BAM/
#
# RESOURCES (per job):
#   - Full run: 32 CPUs, 256GB RAM, 150 hours (6.25 days) - ALL samples
#   - Test run: 8 CPUs, 32GB RAM, 12 hours - Single sample only
#   - All samples processed in a single massive job (normal mode)
#
# NOTES:
#   - Dry-run mode: Use --dry-run to preview without submitting
#   - Normal mode: Processes ALL samples in the dataset in one job
#   - Test mode: Use --sample <NAME> to test with a single sample
#   - Results are automatically copied back to scratch storage
#   - Plots are generated using R (ggplot2, dplyr, gridExtra)
#   - Coverage analysis uses fast mode with 250bp bins (faster and smaller output)
#
# AUTHOR: Generated with AI assistance
# DATE: 2025-01-05
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

main() {
    # Check for command line flags
    local dry_run=false
    local test_sample=""
    
    while [[ $# -gt 0 ]]; do
        case $1 in
            --dry-run)
                dry_run=true
                echo "🧪 DRY RUN MODE - No jobs will be submitted"
                echo ""
                shift
                ;;
            --sample)
                test_sample="$2"
                echo "🧪 TEST MODE - Processing single sample: $test_sample"
                echo ""
                shift 2
                ;;
            *)
                echo "Unknown option: $1"
                echo "Usage: bash BAM_mosdepth_pipeline.sh [--dry-run] [--sample <SAMPLE_NAME>]"
                exit 1
                ;;
        esac
    done
    
    show_welcome
    
    local dataset_name=""
    
    # Keep trying until valid dataset is found or user cancels
    while [ -z "$dataset_name" ]; do
        local user_input=$(get_dataset_name)
        
        # Find actual dataset folder (case-insensitive)
        dataset_name=$(find_dataset_folder "$user_input")
        
        if [ -z "$dataset_name" ]; then
            echo ""
            echo "❌ Dataset '$user_input' not found"
            echo ""
            
            # Let user choose from available datasets
            dataset_name=$(choose_from_available_datasets)
            
            # If user chose to retry (empty string returned), continue the loop
            if [ -z "$dataset_name" ]; then
                continue  # Go back to prompting for dataset name
            fi
        fi
    done
    
    local bam_path="/QRISdata/Q8367/WGS_Reference_Panel/${dataset_name}/4.BAM"
    
    validate_path "$bam_path" "$dataset_name"
    
    # Confirmation checkpoint
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "Dataset: $dataset_name"
    echo "Path: $bam_path"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo ""
    
    if ! confirm "Proceed with this dataset?"; then
        echo "✅ Cancelled by user."
        exit 0
    fi
    
    echo ""
    create_sample_list_and_submit "$dataset_name" "$bam_path" "$dry_run" "$test_sample"
}

show_welcome() {
    clear
    echo ""
    echo "╔═══════════════════════════════════════════════════════════╗"
    echo "║       📊 BAM MOSDEPTH ANALYSIS PIPELINE 📊               ║"
    echo "╚═══════════════════════════════════════════════════════════╝"
    echo ""
}

get_dataset_name() {
    local dataset_name
    while true; do
        read -p "Please state your dataset: " dataset_name
        if [ -n "$dataset_name" ]; then
            echo "$dataset_name"
            return
        fi
        echo "❌ Dataset name cannot be empty."
    done
}

find_dataset_folder() {
    local user_input="$1"
    local base_path="/QRISdata/Q8367/WGS_Reference_Panel"
    
    # Try exact match first (case-sensitive)
    if [ -d "$base_path/$user_input" ]; then
        echo "$user_input"
        return
    fi
    
    # Case-insensitive search: iterate through all directories
    if [ -d "$base_path" ]; then
        for dir in "$base_path"/*; do
            if [ -d "$dir" ]; then
                local dirname=$(basename "$dir")
                # Case-insensitive comparison
                if [[ "$dirname" == "$user_input" ]] || [[ "${dirname,,}" == "${user_input,,}" ]]; then
                    echo "$dirname"
                    return
                fi
            fi
        done
    fi
    
    # Not found
    echo ""
}

choose_from_available_datasets() {
    local base_path="/QRISdata/Q8367/WGS_Reference_Panel"
    
    # Check if base path exists
    if [ ! -d "$base_path" ]; then
        echo "❌ Base path does not exist: $base_path" >&2
        echo "   Please check if you're on the correct system." >&2
        echo "" >&2
        echo ""
        return 0
    fi
    
    # Get list of available datasets
    local datasets=()
    while IFS= read -r dir; do
        datasets+=("$dir")
    done < <(find "$base_path" -maxdepth 1 -type d ! -path "$base_path" -exec basename {} \; 2>/dev/null | sort)
    
    if [ ${#datasets[@]} -eq 0 ]; then
        echo "❌ No datasets found in $base_path" >&2
        echo "   Available directories in base path:" >&2
        ls -la "$base_path" >&2
        echo "" >&2
        echo ""
        return 0
    fi
    
    echo "Available datasets:" >&2
    echo "" >&2
    for i in "${!datasets[@]}"; do
        printf "  %2d) %s\n" $((i+1)) "${datasets[$i]}" >&2
    done
    echo "" >&2
    echo "  r) Retry typing dataset name" >&2
    echo "  0) Cancel and exit" >&2
    echo "" >&2
    
    while true; do
        read -p "Select a number: " choice
        
        # Retry option - return empty to continue loop in main
        if [ "$choice" = "r" ] || [ -z "$choice" ]; then
            echo ""  # Return empty string to main loop (only stdout, no messages)
            return 0
        fi
        
        # Cancel option - exit completely
        if [ "$choice" = "0" ]; then
            echo ""
            echo "✅ Exiting..." >&2
            exit 0
        fi
        
        # Validate number
        if [[ "$choice" =~ ^[0-9]+$ ]] && [ "$choice" -ge 1 ] && [ "$choice" -le ${#datasets[@]} ]; then
            local selected="${datasets[$((choice-1))]}"
            echo "" >&2
            echo "✓ Selected: $selected" >&2
            echo "" >&2
            echo "$selected"  # Return ONLY the selected dataset name (last echo to stdout)
            return 0
        else
            echo "❌ Invalid selection. Please choose 1-${#datasets[@]}, 'r' to retry, or '0' to cancel" >&2
        fi
    done
}

validate_path() {
    local path="$1"
    local dataset_name="$2"
    
    echo "Validating: $path"
    
    if [ ! -d "$path" ]; then
        echo ""
        echo "❌ ERROR: Dataset not found or not in the correct location"
        echo ""
        echo "Expected path: $path"
        echo "Dataset: $dataset_name"
        echo ""
        echo "Please check:"
        echo "  1. Dataset name is correct (case-insensitive)"
        echo "  2. Path exists on the system"
        echo "  3. You have access to /QRISdata/Q8367/WGS_Reference_Panel/"
        exit 1
    fi
    
    echo "✓ Path exists"
    echo ""
}

create_sample_list_and_submit() {
    local dataset_name="$1"
    local bam_path="$2"
    local dry_run="$3"
    local test_sample="$4"
    
    # STEP 1: Scan 4.BAM folder and count samples
    # This list will be used by the SLURM script to process all BAM files
    echo "Scanning for BAM files in $bam_path..."
    local sample_list_file="${SCRIPT_DIR}/bam_list_${dataset_name}_$(date +%Y%m%d_%H%M%S).txt"
    local sample_count=0
    
    > "$sample_list_file"  # Create empty file to start
    
    # Check if test mode (single sample)
    if [ -n "$test_sample" ]; then
        echo "🧪 Test mode: Looking for sample: $test_sample"
        
        # Try to find the sample directory (case-insensitive)
        local found_sample_dir=""
        for sample_dir in "$bam_path"/*; do
            if [ -d "$sample_dir" ]; then
                local dirname=$(basename "$sample_dir")
                # Case-insensitive comparison
                if [[ "${dirname,,}" == "${test_sample,,}" ]]; then
                    found_sample_dir="$sample_dir"
                    break
                fi
            fi
        done
        
        if [ -z "$found_sample_dir" ]; then
            echo "❌ Sample '$test_sample' not found in $bam_path"
            echo ""
            echo "Available samples (first 10):"
            ls -1 "$bam_path" | head -n 10
            exit 1
        fi
        
        # Get recalibrated BAM file from the found sample directory
        local bam_file=$(find "$found_sample_dir" -maxdepth 1 -name "*_recal.bam" -type f | head -n 1)
        if [ -n "$bam_file" ]; then
            echo "$bam_file" >> "$sample_list_file"
            sample_count=1
            echo "✓ Found sample: $test_sample"
            echo "  Recal BAM: $bam_file"
        else
            echo "❌ No recalibrated BAM file (*_recal.bam) found in sample directory: $found_sample_dir"
            exit 1
        fi
    else
        # Normal mode: find all recalibrated BAM files in subdirectories (one per sample directory)
        for sample_dir in "$bam_path"/*; do
            if [ -d "$sample_dir" ]; then
                # Get the recalibrated BAM file from this sample directory
                local bam_file=$(find "$sample_dir" -maxdepth 1 -name "*_recal.bam" -type f | head -n 1)
                if [ -n "$bam_file" ]; then
                    echo "$bam_file" >> "$sample_list_file"  # Write full path to list
                    sample_count=$((sample_count + 1))
                fi
            fi
        done
        
        if [ $sample_count -eq 0 ]; then
            echo "❌ No recalibrated BAM files (*_recal.bam) found"
            rm -f "$sample_list_file"
            exit 1
        fi
    fi
    
    local array_max=$((sample_count - 1))
    
    if [ -n "$test_sample" ]; then
        echo "✓ Found 1 sample for testing"
        echo "✓ Sample: $test_sample"
    else
        echo "✓ Found $sample_count samples in 4.BAM folder"
    fi
    echo ""
    head -n 3 "$sample_list_file" | nl
    echo ""
    
    # Determine message based on mode
    local job_type="samples"
    if [ -n "$test_sample" ]; then
        job_type="TEST sample"
    fi
    
    if [ "$dry_run" = true ]; then
        echo ""
        echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
        echo "🧪 DRY RUN - Would submit job but skipping..."
        echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
        echo ""
        echo "Would submit with:"
        echo "  - Dataset: $dataset_name"
        echo "  - Samples: $sample_count"
        if [ -n "$test_sample" ]; then
            echo "  - Test sample: $test_sample"
        fi
        echo "  - BAM list: $sample_list_file"
        echo ""
        submit_slurm_job "$dataset_name" "$sample_list_file" "$sample_count" "$dry_run" "$test_sample"
        echo ""
        echo "✅ Dry run complete. Check the generated SLURM script for details."
        echo "Sample list: $sample_list_file"
    elif confirm "Submit SLURM job for $sample_count $job_type?"; then
        submit_slurm_job "$dataset_name" "$sample_list_file" "$sample_count" "$dry_run" "$test_sample"
    else
        echo "✅ Cancelled. Sample list: $sample_list_file"
        exit 0
    fi
}

submit_slurm_job() {
    local dataset_name="$1"
    local sample_list_file="$2"
    local sample_count="$3"
    local dry_run="$4"
    local test_sample="$5"
    
    # Generate unique SLURM script filename with timestamp
    local slurm_script="${SCRIPT_DIR}/mosdepth_${dataset_name}_$(date +%Y%m%d_%H%M%S).sh"
    
    # Adjust resources based on mode (test vs full run)
    # Full run: 36 CPUs for parallel mosdepth (6 samples x 6 threads)
    local cpus=36
    local mem="256G"
    local time="150:00:00"
    local parallel_jobs=6
    
    if [ -n "$test_sample" ]; then
        # Much lighter resources for single sample test
        cpus=8
        mem="64G"
        time="12:00:00"
        parallel_jobs=1  # Only 1 sample in test mode
        echo "🧪 TEST MODE: Using reduced resources: $cpus CPUs, $mem RAM, $time time"
    else
        echo "Full run: Using 36 CPUs for parallel mosdepth (6 samples in parallel x 6 threads each)"
    fi
    
    # STEP 2: Generate the SLURM script that will run on compute node
    # This heredoc uses single quotes to prevent variable expansion (we replace them later with sed)
    cat > "$slurm_script" << 'EOF'
#!/bin/bash
#SBATCH --job-name=mosdepth_merge
#SBATCH --output=mosdepth_merge_%j.out
#SBATCH --error=mosdepth_merge_%j.err
#SBATCH --cpus-per-task=CPUS_VALUE
#SBATCH --mem=MEM_VALUE
#SBATCH --time=TIME_VALUE
#SBATCH --account=a_qaafi_cas
#SBATCH --partition=general
#SBATCH -o /scratch/user/uqpha1/logs/mosdepth/mosdepth_merge_%j.output
#SBATCH -e /scratch/user/uqpha1/logs/mosdepth/mosdepth_merge_%j.error


# ============================================================================
# SLURM SCRIPT: Mosdepth Coverage Analysis and Merging
# ============================================================================
# This script runs on the compute node and processes ALL samples together
# BAM_LIST and SAMPLE_COUNT are replaced by sed before submission
# ============================================================================

# Variables set by the wrapper script (replaced via sed):
BAM_LIST="BAM_LIST_FILE"   # Full path to file containing all BAM file paths
SAMPLE_COUNT="SAMPLE_COUNT" # Number of samples to process
PARALLEL_JOBS="PARALLEL_JOBS_VALUE" # Number of parallel jobs (6 for full run, 1 for test)

echo "=========================================="
echo "Processing $SAMPLE_COUNT samples"
echo "BAM list file: $BAM_LIST"
echo "Parallel jobs: $PARALLEL_JOBS"
echo "=========================================="
echo ""

# Verify the BAM list file exists
echo "Verifying BAM list file exists..."
if [ ! -f "$BAM_LIST" ]; then
    echo "ERROR: BAM list file not found: $BAM_LIST"
    exit 1
fi
echo "✓ BAM list file found"
echo ""
echo "First 3 BAM files in list:"
head -n 3 "$BAM_LIST" | nl
echo "=========================================="
echo ""

# Setup environment: Load conda and R plotting environment
module load miniforge
source $ROOTMINIFORGE/etc/profile.d/conda.sh
conda activate rplot

# Change to TMPDIR (local SSD) for fast I/O
cd $TMPDIR

# ============================================================================
# STEP 1: Copy all BAM files to TMPDIR (local SSD for fast processing)
# ============================================================================
# Copy each BAM file and its index from the list to TMPDIR
echo "Copying all BAM files to TMPDIR..."
echo "Reading from: $BAM_LIST"
line_count=0
while IFS= read -r bam_file; do
    if [ -n "$bam_file" ]; then
        line_count=$((line_count + 1))
        echo "[$line_count] Copying: $bam_file"
        rsync -rhivPt "$bam_file" "$TMPDIR"/
        
        # Copy BAI index if it exists (required by mosdepth)
        bai_file="${bam_file}.bai"
        if [ -f "$bai_file" ]; then
            rsync -rhivPt "$bai_file" "$TMPDIR"/
        fi
    fi
done < "$BAM_LIST"

echo ""
echo "✓ Copied $line_count BAM files to TMPDIR"
echo ""

# ============================================================================
# STEP 2: Run mosdepth on each sample in parallel
# ============================================================================
# mosdepth generates regional coverage files in 250bp bins
# Output: .regions.bed.gz files (not per-base due to -n and --fast-mode flags)
echo "Running mosdepth for all samples in parallel (6 samples at a time)..."

# Track successful and failed samples
successful_samples=0
failed_samples=0
start_time_total=$(date +%s)

# Define a function to process one sample
process_sample() {
    local bam_file="$1"
    bam_filename=$(basename "$bam_file")
    # Extract base sample name (remove _recal.bam extension)
    sample_name="${bam_filename%_recal.bam}"
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting: $sample_name"
    start_time=$(date +%s)
    
    # Run mosdepth with error checking
    # -n: no per-nucleotide output
    # --fast-mode: skip coverage per-nucleotide and depth-related metrics
    # --by 250: quantize coverage to every 250th base
    # -t 6: Use 6 threads per sample
    if mosdepth -n --fast-mode --by 250 -t 6 "$sample_name" "$bam_filename"; then
        end_time=$(date +%s)
        duration=$((end_time - start_time))
        minutes=$((duration / 60))
        seconds=$((duration % 60))
        
        # Verify output files were created
        if [ -f "${sample_name}.regions.bed.gz" ]; then
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] ✓ SUCCESS: $sample_name completed in ${minutes}m ${seconds}s"
            echo 1  # Return success indicator
        else
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] ❌ ERROR: Output file not found for $sample_name"
            echo 0  # Return failure indicator
        fi
    else
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] ❌ ERROR: mosdepth failed for $sample_name"
        echo 0  # Return failure indicator
    fi
}

# Export the function so xargs can use it
export -f process_sample

# Run mosdepth in parallel using xargs
# -P $PARALLEL_JOBS: Process samples in parallel (6 for full run, 1 for test)
# -I {}: Replace {} with the line from input
# Capture exit codes
echo "Running with $PARALLEL_JOBS parallel jobs (using $((PARALLEL_JOBS * 6)) CPUs)"
results=$(while IFS= read -r bam_file; do
    if [ -n "$bam_file" ]; then
        echo "$bam_file"
    fi
done < "$BAM_LIST" | xargs -n 1 -P $PARALLEL_JOBS -I {} bash -c 'process_sample "{}"')

# Count successes and failures from the output
# Use wc -l to count lines containing SUCCESS/ERROR
successful_samples=$(echo "$results" | grep "SUCCESS" | wc -l)
failed_samples=$(echo "$results" | grep "ERROR" | wc -l)
# Ensure we have valid numbers (handle empty case)
[ -z "$successful_samples" ] && successful_samples=0
[ -z "$failed_samples" ] && failed_samples=0

# Print summary
end_time_total=$(date +%s)
duration_total=$((end_time_total - start_time_total))
minutes_total=$((duration_total / 60))
seconds_total=$((duration_total % 60))

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Mosdepth Summary"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Total samples processed: $((successful_samples + failed_samples))"
echo "✓ Successful: $successful_samples"
echo "❌ Failed: $failed_samples"
echo "Total time: ${minutes_total}m ${seconds_total}s"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

# Exit if any samples failed
if [ $failed_samples -gt 0 ]; then
    echo "ERROR: Some samples failed during mosdepth analysis"
    exit 1
fi

# ============================================================================
# STEP 3: Merge all regions coverage files
# ============================================================================
# Combine coverage from all samples by averaging at each position
# Note: With -n --fast-mode --by 250, mosdepth generates .regions.bed.gz files
echo "Merging mosdepth results..."
echo "Merging regions coverage files (250bp bins)..."
zcat *.regions.bed.gz | sort -k1,1 -k2,2n | \
awk 'BEGIN {chr=""; start=0; end=0; total_coverage=0; count=0}
     {
       if (chr != $1 || start != $2) {
         if (chr != "") {
           print chr, start, end, total_coverage/count
         }
         chr=$1; start=$2; end=$3; total_coverage=$4; count=1
       } else {
         total_coverage += $4; count++
       }
     }
     END {if (chr != "") print chr, start, end, total_coverage/count}' | \
gzip > merged_coverage.bed.gz

echo ""

# ============================================================================
# STEP 4: Generate visualization plots with R
# ============================================================================
# Creates two types of plots:
#  1. Faceted plot: All chromosomes in one grid view
#  2. Individual plot: Each chromosome in its own subplot
echo "Generating merged coverage plot..."
echo "Checking if merged_coverage.bed.gz exists and has data..."

# Check if merged file exists
if [ ! -f "merged_coverage.bed.gz" ]; then
    echo "ERROR: merged_coverage.bed.gz not found"
    exit 1
fi

# Check if merged file has data
line_count=$(zcat merged_coverage.bed.gz | wc -l)
echo "✓ Merged coverage file found with $line_count lines"
echo ""

# Export SAMPLE_COUNT for R to read (value will be replaced by sed)
export SAMPLE_COUNT_VAR="SAMPLE_COUNT_VALUE"

echo "Starting R script for plotting..."
# Generate both faceted and individual chromosome plots using ggplot2
Rscript -e "
library(ggplot2)
library(dplyr)

# Read merged coverage file
coverage_data <- read.table('merged_coverage.bed.gz', header=FALSE, 
                           col.names=c('chr', 'start', 'end', 'coverage'))

# Filter for chromosomes Chr00-Chr17 (18 chromosomes)
chr_list <- paste0('Chr', sprintf('%02d', 0:17))
coverage_data <- coverage_data[coverage_data\$chr %in% chr_list, ]

# Convert chromosome to factor for proper ordering
coverage_data\$chr <- factor(coverage_data\$chr, levels=chr_list)

# Convert position to MegaBasepairs
coverage_data\$pos_mb <- coverage_data\$start / 1000000

# Get sample count from bash variable
sample_count <- Sys.getenv('SAMPLE_COUNT_VAR')

# 1. Create faceted plot (all chromosomes in one plot)
p_faceted <- ggplot(coverage_data, aes(x=pos_mb, y=coverage)) +
  geom_line(color='steelblue', alpha=0.7, linewidth=0.5) +
  facet_wrap(~chr, scales='free_x', ncol=3) +
  theme_minimal() +
  theme(strip.text = element_text(size=10),
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8)) +
  labs(title='Average Coverage Across Chromosomes - DATASET_NAME',
       x='Position (Mbp)',
       y='Average Coverage (X)',
       subtitle=paste('Merged coverage from', sample_count, 'samples'))

ggsave('merged_coverage_faceted.png', p_faceted, width=15, height=12, dpi=300)
cat('Faceted plot saved: merged_coverage_faceted.png\n')

# 2. Create individual line plots for each chromosome in one combined plot
library(gridExtra)

# Create list to store individual plots
plot_list <- list()

for (chr in chr_list) {
  chr_data <- coverage_data[coverage_data\$chr == chr, ]
  
  if (nrow(chr_data) > 0) {
    p <- ggplot(chr_data, aes(x=pos_mb, y=coverage)) +
      geom_line(color='steelblue', alpha=0.7, linewidth=0.5) +
      theme_minimal() +
      theme(axis.text.x = element_text(size=8),
            axis.text.y = element_text(size=8),
            plot.title = element_text(size=10, hjust=0.5),
            plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), 'cm')) +
      labs(title=chr,
           x='Position (Mbp)',
           y='Coverage (X)') +
      ylim(0, 200)
    
    plot_list[[chr]] <- p
  }
}

# Combine all plots into one grid
combined_plot <- do.call(grid.arrange, c(plot_list, ncol=3))
ggsave('individual_chromosomes_combined.png', combined_plot, width=20, height=24, dpi=300)
cat('Combined individual chromosome plots saved: individual_chromosomes_combined.png\n')

cat('All plots completed: 1 faceted + 1 combined individual chromosome plot\n')
"

# Verify R script completed successfully
if [ $? -eq 0 ]; then
    echo "✓ R script completed successfully"
    
    # Check if output files were created
    if [ -f "merged_coverage_faceted.png" ]; then
        echo "✓ Faceted plot created: merged_coverage_faceted.png"
    else
        echo "❌ ERROR: merged_coverage_faceted.png not found"
        exit 1
    fi
    
    if [ -f "individual_chromosomes_combined.png" ]; then
        echo "✓ Individual chromosomes plot created: individual_chromosomes_combined.png"
    else
        echo "❌ ERROR: individual_chromosomes_combined.png not found"
        exit 1
    fi
    echo ""
else
    echo "❌ ERROR: R script failed"
    exit 1
fi

echo ""

# ============================================================================
# STEP 5: Copy results back to scratch storage
# ============================================================================
# Copy all output files from TMPDIR back to permanent scratch storage
echo "Copying results to scratch..."

# Create output directory and copy results back to scratch
OUTPUT_DIR="/scratch/user/uqpha1/mosdepth/${dataset_name}_merged"
mkdir -p "$OUTPUT_DIR"

echo "Copying merged results to $OUTPUT_DIR"
rsync -rhivPt merged_coverage.bed.gz "$OUTPUT_DIR"/
rsync -rhivPt merged_coverage_faceted.png "$OUTPUT_DIR"/
rsync -rhivPt individual_chromosomes_combined.png "$OUTPUT_DIR"/

echo ""
echo "✓ Complete: Merged analysis for $SAMPLE_COUNT samples"
echo "✓ Results saved to: $OUTPUT_DIR"
EOF

    # STEP 3: Replace placeholders with actual values
    # The heredoc above used placeholders to avoid variable expansion
    # Now we replace them with the real values from the main script
    # Note: Using -i with backup for compatibility across Linux/macOS
    
    # First, replace the standalone variable assignments and values
    sed -i.bak "s|CPUS_VALUE|$cpus|g" "$slurm_script"
    sed -i.bak "s|MEM_VALUE|$mem|g" "$slurm_script"
    sed -i.bak "s|TIME_VALUE|$time|g" "$slurm_script"
    sed -i.bak "s|DATASET_NAME|$dataset_name|g" "$slurm_script"
    sed -i.bak "s|PARALLEL_JOBS_VALUE|$parallel_jobs|g" "$slurm_script"
    
    # Replace BAM_LIST_FILE (must be before SAMPLE_COUNT to avoid conflicts)
    sed -i.bak "s|BAM_LIST_FILE|$sample_list_file|g" "$slurm_script"
    
    # Replace SAMPLE_COUNT but ONLY in specific contexts to avoid corrupting other strings
    # Replace in variable assignment: SAMPLE_COUNT="SAMPLE_COUNT" -> SAMPLE_COUNT="50"
    sed -i.bak "s|SAMPLE_COUNT=\"SAMPLE_COUNT\"|SAMPLE_COUNT=\"$sample_count\"|g" "$slurm_script"
    
    # Replace SAMPLE_COUNT_VALUE placeholder
    sed -i.bak "s|SAMPLE_COUNT_VALUE|$sample_count|g" "$slurm_script"
    
    # Replace in echo statements using the specific pattern (with escaped $)
    sed -i.bak "s|Processing \$SAMPLE_COUNT samples|Processing $sample_count samples|g" "$slurm_script"
    sed -i.bak "s|Merged coverage from \$SAMPLE_COUNT samples|Merged coverage from $sample_count samples|g" "$slurm_script"
    sed -i.bak "s|Merged analysis for \$SAMPLE_COUNT samples|Merged analysis for $sample_count samples|g" "$slurm_script"
    
    # Remove backup file created by sed
    rm -f "${slurm_script}.bak"
    
    # Debug: Show what BAM_LIST was set to
    echo "Debug: BAM_LIST variable set to: $(grep '^BAM_LIST=' "$slurm_script" | head -n 1)"
    
    # Make script executable
    chmod +x "$slurm_script"
    
    echo ""
    # STEP 4: Submit the job or preview it (dry-run)
    if [ "$dry_run" = true ]; then
        # Dry-run mode: Show what would happen without submitting
        echo "🧪 DRY RUN - Would execute: sbatch $slurm_script"
        echo ""
        echo "Generated SLURM script: $slurm_script"
        echo ""
        echo "First 50 lines of generated script:"
        head -n 50 "$slurm_script"
        echo ""
        echo "... (rest of script would run as normal) ..."
    else
        # Normal mode: Actually submit the job to SLURM
        echo "Submitting SLURM job..."
        local job_id=$(sbatch "$slurm_script" | awk '{print $NF}')
        
        if [ $? -eq 0 ]; then
            echo ""
            echo "✅ SUCCESS!"
            echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
            echo "Job ID: $job_id"
            if [ -n "$test_sample" ]; then
                echo "Mode: TEST run (single sample)"
                echo "Test sample: $test_sample"
            else
                echo "Samples: $sample_count (single massive job)"
            fi
            echo "Resources: $cpus CPUs, $mem RAM, $time"
            echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
            echo "Monitor: squeue -j $job_id"
            echo "Cancel:  scancel $job_id"
        else
            echo "❌ Failed to submit job"
            exit 1
        fi
    fi
}

confirm() {
    local message="$1"
    while true; do
        read -p "$message [y/n]: " response
        case $response in
            [Yy]* ) return 0;;
            [Nn]* ) return 1;;
            * ) echo "Please answer y or n.";;
        esac
    done
}

main "$@"
