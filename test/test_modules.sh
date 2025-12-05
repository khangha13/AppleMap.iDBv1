#!/bin/bash
# =============================================================================
# MODULE TESTING SCRIPT
# =============================================================================
# Test script to verify that all modules work correctly

echo "🧪 Testing GATK Pipeline Modules"
echo "================================="
echo ""

# Test Step 1A Functions
echo "📋 Testing Step 1A Functions..."
echo "-------------------------------"

# Source Step 1A functions
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="$(dirname "$SCRIPT_DIR")"
source "$PIPELINE_DIR/modules/step1a/lib/functions.sh" 2>/dev/null

if [ $? -eq 0 ]; then
    echo "✅ Step 1A functions sourced successfully"
    
    # Test if key functions are defined
    if declare -f start_step_timer >/dev/null 2>&1; then
        echo "✅ start_step_timer function defined"
    else
        echo "❌ start_step_timer function not found"
    fi
    
    if declare -f run_fastqc >/dev/null 2>&1; then
        echo "✅ run_fastqc function defined"
    else
        echo "❌ run_fastqc function not found"
    fi
    
    if declare -f run_haplotype_caller >/dev/null 2>&1; then
        echo "✅ run_haplotype_caller function defined"
    else
        echo "❌ run_haplotype_caller function not found"
    fi
else
    echo "❌ Failed to source Step 1A functions"
fi

echo ""

# Test Step 1B Functions
echo "📋 Testing Step 1B Functions..."
echo "-------------------------------"

# Source Step 1B functions
source "$PIPELINE_DIR/modules/step1b/lib/functions.sh" 2>/dev/null

if [ $? -eq 0 ]; then
    echo "✅ Step 1B functions sourced successfully"
    
    # Test if key functions are defined
    if declare -f build_sample_map >/dev/null 2>&1; then
        echo "✅ build_sample_map function defined"
    else
        echo "❌ build_sample_map function not found"
    fi
    
    if declare -f run_genomics_db_import >/dev/null 2>&1; then
        echo "✅ run_genomics_db_import function defined"
    else
        echo "❌ run_genomics_db_import function not found"
    fi
    
    if declare -f get_chromosome_list >/dev/null 2>&1; then
        echo "✅ get_chromosome_list function defined"
    else
        echo "❌ get_chromosome_list function not found"
    fi
else
    echo "❌ Failed to source Step 1B functions"
fi

echo ""

# Test Main Entry Point
echo "📋 Testing Main Entry Point..."
echo "-----------------------------"

# Test if bin/gatk_pipeline.sh exists and is executable
if [ -f "$PIPELINE_DIR/bin/gatk_pipeline.sh" ]; then
    echo "✅ bin/gatk_pipeline.sh exists"
    
    if [ -x "$PIPELINE_DIR/bin/gatk_pipeline.sh" ]; then
        echo "✅ bin/gatk_pipeline.sh is executable"
    else
        echo "⚠️  bin/gatk_pipeline.sh is not executable (fixing...)"
        chmod +x "$PIPELINE_DIR/bin/gatk_pipeline.sh"
        echo "✅ bin/gatk_pipeline.sh is now executable"
    fi
else
    echo "❌ bin/gatk_pipeline.sh not found"
fi

echo ""

# Test Module Structure
echo "📋 Testing Module Structure..."
echo "-----------------------------"

# Check if all required directories exist
required_dirs=(
    "modules/step1a"
    "modules/step1b"
    "lib"
    "config"
    "utils"
)

for dir in "${required_dirs[@]}"; do
    if [ -d "$PIPELINE_DIR/$dir" ]; then
        echo "✅ $dir directory exists"
    else
        echo "❌ $dir directory missing"
    fi
done

echo ""

# Test Required Files
echo "📋 Testing Required Files..."
echo "----------------------------"

required_files=(
    "modules/step1a/lib/functions.sh"
    "modules/step1a/bin/run_step1a.sh"
    "modules/step1b/lib/functions.sh"
    "modules/step1b/bin/run_step1b.sh"
    "lib/logging.sh"
    "lib/slurm.sh"
    "lib/validation.sh"
    "lib/pipeline_common.sh"
    "bin/gatk_pipeline.sh"
)

for file in "${required_files[@]}"; do
    if [ -f "$PIPELINE_DIR/$file" ]; then
        echo "✅ $file exists"
    else
        echo "❌ $file missing"
    fi
done

echo ""
echo "📋 Testing Step 1D Scripts..."
echo "----------------------------"

MASTER_VCF_SCRIPT="$PIPELINE_DIR/modules/step1d/templates/master_vcf_analysis.sh"
if bash -n "$MASTER_VCF_SCRIPT"; then
    echo "✅ master_vcf_analysis.sh passes bash -n syntax check"
else
    echo "❌ master_vcf_analysis.sh failed bash -n syntax check"
    exit 1
fi

INTERACTIVE_WRAPPER="$PIPELINE_DIR/wrappers/interactive/step1d_interactive.sh"
if bash -n "$INTERACTIVE_WRAPPER"; then
    echo "✅ step1d_interactive.sh passes bash -n syntax check"
else
    echo "❌ step1d_interactive.sh failed bash -n syntax check"
    exit 1
fi

echo ""
echo "🎉 Module Testing Complete!"
echo "=========================="
