# 🎯 VCF QC Analysis Pipeline - Quick Guide

**Ready to analyze your VCF files on HPC? Start here!**

---

## 📚 Documentation Guide

Read these files in this order:

### 1️⃣ First Time on HPC
**Read:** `HPC_SETUP_GUIDE.md`
- How to upload files
- How to configure
- How to run
- Complete troubleshooting

### 2️⃣ Quick Reference
**Read:** `START_HERE.md`
- Quick overview
- 3-step quick start
- File descriptions

### 3️⃣ Configuration
**Edit:** `vcf_analysis_config.sh`
- Set your VCF location
- Set your HPC account
- Verify settings

### 4️⃣ Deployment
**Use:** `HPC_DEPLOYMENT_CHECKLIST.md`
- Step-by-step checklist
- Validation steps
- Verification procedures

### 5️⃣ Understanding the Pipeline
**Read:** `VCF_Analysis_WORKFLOW.md`
- Visual diagrams
- Data flow
- Execution modes

---

## ⚡ Super Quick Start

Already know what you're doing?

```bash
# 1. Edit config
nano vcf_analysis_config.sh
# Set: VCF_DIR, SLURM_ACCOUNT, LOG_DIR

# 2. Validate
source vcf_analysis_config.sh && validate_config

# 3. Run
bash master_vcf_analysis_array.sh

# 4. Monitor
squeue -u $USER
```

---

## 🆘 Need Help?

- **Setup issues:** → `HPC_SETUP_GUIDE.md`
- **Quick commands:** → `VCF_Analysis_Quick_Reference.md`
- **Understanding flow:** → `VCF_Analysis_WORKFLOW.md`
- **What changed:** → `HPC_MIGRATION_SUMMARY.md`

---

## ✅ What This Pipeline Does

```
Input: 18 VCF files (Chr00-Chr17)
  ↓
Extract depth → Generate plots → Analyze missingness → Combine results
  ↓                ↓                 ↓                     ↓
TSV file     PDF plots        PNG per chr          Combined PNG
```

**Time:** ~20-26 hours with parallel processing
**Resources:** 18-24 CPUs, 64-128GB RAM

---

## 🧮 Optional PCA Stage (Step 1D+)

Need a quick look at cohort structure after QC? Run PCA mode:

- **Interactive wrapper:** `wrappers/interactive/step1d_interactive.sh --PCA [--remove-relatives]`
- **Batch wrapper:** `wrappers/sbatch/step1d_submit.sh [<dataset>] <vcf_dir> --PCA` (defaults to `<vcf_dir>` basename if omitted)
- **Manual:** `bash modules/step1d/templates/master_vcf_analysis.sh --PCA` (set `STEP1D_REMOVE_RELATIVES=true` for KING 0.125 filtering).

Requirements: `bcftools` and `plink2` available on the compute node plus R packages `ggplot2`, `data.table`, `ragg`, `scales`.  
Outputs land under `${WORK_DIR}/pca_analysis/` (configurable via `STEP1D_PCA_DIR`) and include `pca.eigenvec`, `pca.eigenval`, and ready-to-share PNGs (`pca_PC1_PC2.png`, `pca_scree.png`).
Use `STEP1D_PCA_FORCE_CONCAT=true` to always concatenate per-chrom VCFs; merged detection ignores filenames containing `Chr` by default (set `STEP1D_PCA_MERGED_EXCLUDE_CHR=false` to allow).

> The automation loads `miniforge/25.3.0-3`, `bcftools`, and `plink/2.00a3.6-gcc-11.3.0` modules by default. Adjust those module names (or preload your own) if your HPC environment differs.

---

## 🎯 Critical Files

| File | Purpose |
|------|---------|
| `vcf_analysis_config.sh` | **EDIT THIS FIRST** |
| `master_vcf_analysis_array.sh` | **RUN THIS** (parallel) |
| `HPC_SETUP_GUIDE.md` | **READ FOR SETUP** |
| `HPC_DEPLOYMENT_CHECKLIST.md` | **USE FOR DEPLOYMENT** |

---

## 📦 Files in This Directory

```
1c_VCF_QC/
├── README_FIRST.md                     ← You are here!
├── HPC_SETUP_GUIDE.md                  ⭐ Complete setup guide
├── HPC_DEPLOYMENT_CHECKLIST.md         ✅ Deployment checklist
├── START_HERE.md                       🚀 Quick start
├── vcf_analysis_config.sh              ⚙️  EDIT THIS
├── master_vcf_analysis_array.sh        🎯 RUN THIS (parallel)
├── master_vcf_analysis.sh              Single job version
├── extract_site_mean_DP.sh             Standalone TSV extraction
└── VCF_Analysis_*.md                   Additional documentation
```

---

## 🎓 First Time?

1. Read `HPC_SETUP_GUIDE.md` (10 minutes)
2. Upload to HPC (5 minutes)
3. Edit `vcf_analysis_config.sh` (5 minutes)
4. Follow `HPC_DEPLOYMENT_CHECKLIST.md` (15 minutes)
5. Run pipeline (automated, 20-26 hours)

**Total hands-on time: ~35 minutes**

---

## 💡 Pro Tip

Test with 2 chromosomes first:
```bash
export CHR_START=0
export CHR_END=1
bash master_vcf_analysis.sh
```

Then run full 18 chromosomes:
```bash
bash master_vcf_analysis_array.sh
```

---

**Questions?** See the documentation files above or check the troubleshooting sections.

**Ready?** Start with `HPC_SETUP_GUIDE.md`!

🧬📊✨ Happy analyzing!
