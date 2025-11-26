# Pipeline Deployment Readiness Report

**Date**: 2025-01-11  
**Pipeline**: GATK Pipeline KH v1  
**Status**: ✅ **READY FOR DEPLOYMENT** (with configuration required)

---

## ✅ Deployment Checklist

### **1. Core Components**
- ✅ Main entry point exists: `bin/gatk_pipeline.sh` (executable)
- ✅ All 4 step modules present: step1a, step1b, step1c, step1d
- ✅ All interactive wrappers present: step1a-1d_interactive.sh
- ✅ All SLURM wrappers present: step1a-1d_submit.sh
- ✅ Configuration template exists: `config/environment.template.sh`
- ✅ Common libraries present: logging.sh, validation.sh, slurm.sh, pipeline_common.sh

### **2. Hardcoded Paths** ✅ FIXED
- ✅ Removed hardcoded `/scratch/user/uqpha1/logs/` → Uses `${LOG_BASE_PATH}`
- ✅ Removed hardcoded `/QRISdata/Q8367/WGS_Reference_panel/` → Uses `${RDM_BASE_PATH}`
- ✅ All paths now use configurable variables from `pipeline_config.sh`
- ✅ Repository location is configurable via `PIPELINE_ROOT` so the codebase can live under `$HOME`, RDM, or scratch

### **3. Directory Structure Consistency** ✅ COMPLETE
- ✅ All scripts updated to match new RDM structure:
  - `6.genomicsdb_output` (Step 1B workspace)
  - `7.Consolidated_VCF` (Step 1B output, Step 1C input)
  - `8.Imputated_VCF_BEAGLE` (Step 1C output)
  - `9.Imputation_QC` (Step 1D output)
  - `9.Metrics` (Step 1A per-sample metrics)

### **4. Configuration System**
- ✅ Environment template provided: `config/environment.template.sh`
- ✅ Pipeline config uses environment variables with defaults
- ✅ All paths configurable via environment file

---

## ⚙️ Pre-Deployment Configuration Required

### **CRITICAL: User must create `config/environment.sh`**

Copy the template and configure:

```bash
cp config/environment.template.sh config/environment.sh
nano config/environment.sh
```

### **Required Settings** (⚠️ Must be configured):

1. **Storage Paths**:
   ```bash
   PIPELINE_ROOT="$HOME/GATK_Pipeline_KH_v1"          # Or absolute path on RDM/scratch
   PIPELINE_RDM_BASE="/QRISdata/QXXXX/PROJECT_NAME"   # Update QXXXX and PROJECT_NAME
   PIPELINE_SCRATCH_BASE="/scratch/user/USERNAME"     # Update USERNAME
   PIPELINE_LOG_BASE="/scratch/user/USERNAME/logs"    # Update USERNAME or point to RDM
   ```

2. **Reference Data**:
   ```bash
   PIPELINE_REFERENCE_DIR="/path/to/reference_directory"
   PIPELINE_REFERENCE_FASTA="${PIPELINE_REFERENCE_DIR}/reference.fasta"
   PIPELINE_KNOWN_SITES_VCF="${PIPELINE_REFERENCE_DIR}/known_sites.vcf.gz"
   PIPELINE_ADAPTER_FASTA="/path/to/adapters/TruSeq3-PE.fa"
   ```

3. **SLURM Configuration**:
   ```bash
   PIPELINE_SLURM_ACCOUNT="your_account_name"          # Update account
   PIPELINE_SLURM_PARTITION="general"                 # Verify partition name
   ```

### **Optional Settings** (Review defaults):
- Step-specific CPU/memory/time limits
- Array job limits
- Verbose logging
- Dry-run mode

---

## 🚀 Deployment Steps

### **1. Upload to HPC**
```bash
rsync -avz GATK_Pipeline_KH_v1/ user@hpc:/path/to/project/GATK_Pipeline_KH_v1/
```

### **2. Configure Environment**
```bash
ssh user@hpc
cd /path/to/project/GATK_Pipeline_KH_v1
cp config/environment.template.sh config/environment.sh
nano config/environment.sh  # Edit required settings
```

### **3. Verify Configuration**
```bash
# Test that config loads correctly
source config/pipeline_config.sh
echo "RDM Base: ${RDM_BASE_PATH}"
echo "Scratch Base: ${SCRATCH_BASE_PATH}"
echo "Log Base: ${LOG_BASE_PATH}"
```

### **4. Test Entry Point**
```bash
# Test help
bin/gatk_pipeline.sh --help

# Test dry-run (if dataset exists)
bin/gatk_pipeline.sh -d test_dataset --dry-run
```

### **5. Verify Run-Specific Log Folder**
- Every submission creates `${LOG_BASE_PATH}/${DATASET}_<YYYYMMDD>/`.
- Slurm stdout/stderr and the structured `logging.sh` file for `gatk_pipeline.sh` are streamed into this folder automatically.
- When using `--submit`/interactive mode, the self-submission step points Slurm’s `-o/-e` flags at the same folder, so all master logs are co-located regardless of whether the repo lives under `$HOME`, RDM, or scratch.
- Interactive runs now stop to ask whether pending steps (1B/1C/1D) should be submitted; if you decline, no jobs are launched. Completed datasets exit immediately with a friendly confirmation.

---

## 📋 Known Limitations

1. **Documentation Files**: Some docs still contain example paths (`/QRISdata/Q8367/`, `/scratch/user/uqpha1/`)
   - **Impact**: Low - these are documentation/examples only
   - **Action**: Can be updated later, doesn't affect functionality

2. **Archive Files**: Old scripts in `archive/` contain hardcoded paths
   - **Impact**: None - these are archived, not used
   - **Action**: No action needed

3. **Utility Scripts**: `utils/mosdepth/` may contain hardcoded paths
   - **Impact**: Low - utility script, not core pipeline
   - **Action**: Can be updated separately if needed

---

## ✅ Deployment Readiness Summary

| Component | Status | Notes |
|-----------|--------|-------|
| Core Scripts | ✅ Ready | All paths configurable |
| Configuration | ✅ Ready | Template provided |
| Directory Structure | ✅ Ready | Consistent across all scripts |
| Hardcoded Paths | ✅ Fixed | All removed from active scripts |
| Documentation | ⚠️ Partial | Some example paths remain (non-critical) |
| **Overall** | ✅ **DEPLOYABLE** | Requires `config/environment.sh` setup |

---

## 🎯 Post-Deployment Validation

After deployment, verify:

1. **Configuration loads**:
   ```bash
   source config/pipeline_config.sh && echo "Config OK"
   ```

2. **Paths are correct**:
   ```bash
   echo "RDM: ${RDM_BASE_PATH}"
   echo "Scratch: ${SCRATCH_BASE_PATH}"
   echo "Logs: ${LOG_BASE_PATH}"
   ```

3. **Entry point works**:
   ```bash
   bin/gatk_pipeline.sh --help
   ```

4. **Interactive mode works**:
   ```bash
   bin/gatk_pipeline.sh  # Should enter interactive mode
   ```

---

## 📝 Notes

- The pipeline is **fully deployable** after creating `config/environment.sh`
- All active scripts use configurable paths
- Documentation examples can be updated post-deployment
- Archive files are intentionally left as-is (historical reference)
- All automation must call `bin/gatk_pipeline.sh`; the legacy entry point without the `.sh` suffix has been removed.

**Recommendation**: ✅ **APPROVED FOR DEPLOYMENT**

