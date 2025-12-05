# VCF Analysis Pipeline - Quick Reference

## 🚀 Quick Commands

### Full Pipeline (One Command)

```bash
# Interactive mode
cd /path/to/Scripts
bash master_vcf_analysis.sh

# SLURM array mode (recommended)
bash master_vcf_analysis_array.sh
```

### Configuration

```bash
# Edit configuration
nano vcf_analysis_config.sh

# View current config
source vcf_analysis_config.sh && show_config

# Validate config
source vcf_analysis_config.sh && validate_config
```

### Monitor SLURM Jobs

```bash
# View all your jobs
squeue -u $USER

# View specific job
squeue -j JOB_ID

# View array task (e.g., Chr05 is task 5)
squeue -j JOB_ID_5

# Cancel job
scancel JOB_ID

# Check job efficiency after completion
seff JOB_ID
```

### Check Outputs

```bash
# List all outputs
ls -lh SNP_site_meanDP.tsv
ls -lh depth_pdfs/*.pdf | wc -l
ls -lh Chr*_missingness_vs_depth.png | wc -l
ls -lh all_chromosomes_missingness_vs_depth.png

# Check TSV file
head SNP_site_meanDP.tsv
wc -l SNP_site_meanDP.tsv

# Check which chromosomes are complete
for i in {0..17}; do
  chr=$(printf "Chr%02d" $i)
  if [ -f "${chr}_missingness_vs_depth.png" ]; then
    echo "✓ $chr complete"
  else
    echo "✗ $chr missing"
  fi
done
```

---

## 📋 Individual Components

### Extract TSV Only

```bash
# Standalone execution
bash extract_site_mean_DP.sh

# Output: SNP_site_meanDP.tsv
```

### Depth Plots Only

```bash
# Requires: SNP_site_meanDP.tsv
Rscript Average_depth_per_Chromosome.R

# Output: depth_pdfs/Chr*.pdf
```

### Missingness for One Chromosome

```bash
# Single chromosome
Rscript per_site_missingness_vs_depth.R Chr05.vcf.gz Chr05_missingness.png

# Batch processing
for i in {0..17}; do
  chr=$(printf "Chr%02d" $i)
  Rscript per_site_missingness_vs_depth.R ${chr}.vcf.gz ${chr}_missingness.png
done
```

### Combine Plots

```bash
# Requires: Chr*_missingness_vs_depth.png
Rscript combine_per_site_missingness_plots.R

# Output: all_chromosomes_missingness_vs_depth.png
```

---

## 🔧 Common Customizations

### Change VCF File Pattern

```bash
# In vcf_analysis_config.sh
export VCF_PATTERN="MyData_Chr%02d_filtered.vcf.gz"
```

### Change Chromosome Range

```bash
# In vcf_analysis_config.sh
export CHR_START=1    # Start from Chr01
export CHR_END=12     # End at Chr12
```

### Increase Memory for Large Files

```bash
# In vcf_analysis_config.sh
export MISSINGNESS_RESOURCES="8|128G|24:00:00"  # 8 CPUs, 128GB RAM, 24 hours
```

### Change Output Directory

```bash
# In vcf_analysis_config.sh
export WORK_DIR="/scratch/user/username/analysis_output"
```

---

## 🐛 Quick Troubleshooting

| Problem | Quick Fix |
|---------|-----------|
| VCF not found | Check `VCF_DIR` and `VCF_PATTERN` in config |
| bcftools not found | `module load bcftools` or `conda install bcftools` |
| R package missing | `Rscript -e "install.packages('PACKAGE_NAME')"` |
| Out of memory | Increase RAM in `*_RESOURCES` variables |
| Permission denied | `chmod +x *.sh` |
| Job pending too long | Check cluster load with `sinfo` |

---

## 📊 Expected Outputs

```
Working Directory/
├── SNP_site_meanDP.tsv                          # ~100MB-2GB (depends on genome size)
├── depth_pdfs/                                  # 18 PDFs
│   ├── Chr00_depth.pdf
│   ├── Chr01_depth.pdf
│   └── ... (through Chr17)
├── Chr00_missingness_vs_depth.png              # 18 PNGs
├── Chr01_missingness_vs_depth.png
├── ... (through Chr17)
└── all_chromosomes_missingness_vs_depth.png    # Combined plot
```

---

## ⏱️ Time Estimates

| Step | Small Dataset | Large Dataset |
|------|--------------|---------------|
| TSV Extraction | 1-4 hours | 12-24 hours |
| Depth Plots | 5-30 min | 1-4 hours |
| Missingness (per chr) | 1-3 hours | 6-12 hours |
| Combine Plots | 2-5 min | 5-15 min |

**Small:** <50 samples, <5M variants per chromosome  
**Large:** >100 samples, >10M variants per chromosome

---

## 📁 File Sizes (Approximate)

| File | Small | Large |
|------|-------|-------|
| SNP_site_meanDP.tsv | 50-200 MB | 500MB-2GB |
| Chr*_depth.pdf | 500KB-2MB each | 2-5MB each |
| Chr*_missingness.png | 2-5MB each | 5-10MB each |
| all_chromosomes.png | 15-30MB | 30-60MB |

---

## 🔄 Re-run Failed Jobs

```bash
# Find failed chromosome
for i in {0..17}; do
  chr=$(printf "Chr%02d" $i)
  [ ! -f "${chr}_missingness_vs_depth.png" ] && echo $i
done

# Re-run specific chromosomes (e.g., 3, 7, 15)
sbatch --array=3,7,15 job_missingness_array.sh
```

---

## 💡 Pro Tips

1. **Start Small**: Test with 1-2 chromosomes first
   ```bash
   export CHR_START=0
   export CHR_END=1
   ```

2. **Check Resources**: Monitor first job to adjust resources
   ```bash
   seff JOB_ID
   ```

3. **Use Array Jobs**: Much faster for missingness plots
   ```bash
   bash master_vcf_analysis_array.sh  # Not master_vcf_analysis.sh
   ```

4. **Keep Logs**: Don't delete log files until analysis is complete
   ```bash
   ls -lh logs/vcf_analysis/
   ```

5. **Backup Outputs**: TSV extraction can take a long time
   ```bash
   cp SNP_site_meanDP.tsv SNP_site_meanDP.tsv.backup
   ```

---

## 📞 Getting Help

```bash
# Check logs
tail -50 logs/vcf_analysis/missingness_JOBID_5.err

# Validate environment
which bcftools
which Rscript
conda list | grep -E "vcfR|ggplot2|data.table"

# Test R script
Rscript --version
Rscript -e "library(vcfR); library(ggplot2)"
```

