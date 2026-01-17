---
name: SNPs-per-sample count table
overview: Add a Step1D-style HPC-friendly utility that scans a directory for chr VCFs, counts per-sample called biallelic PASS SNP genotypes per chromosome, and writes a single R-friendly table containing both per-chromosome rows and per-sample TOTAL rows.
todos:
  - id: inspect-existing-utils
    content: Review existing `utils/` scripts for conventions (SBATCH headers, argument parsing, output naming) to match style.
    status: pending
  - id: design-io-contract
    content: "Define CLI + env var overrides: VCF_DIR/WORK_DIR/OUT_FILE, detection rule (*chr*), dry-run, module name, and filter flags (PASS + biallelic SNPs)."
    status: pending
  - id: implement-counting
    content: "Implement per-VCF loop: filter to PASS biallelic SNPs, run `bcftools stats -s -`, parse PSC, compute n_called, emit per-chr rows."
    status: pending
  - id: aggregate-totals
    content: Aggregate per-sample totals across chr rows and append TOTAL rows to the same output file.
    status: pending
  - id: hpc-hardening
    content: Add PIPELINE_ROOT-safe sourcing, `module purge`, exact module load string, and robust error handling/logging aligned with Step1d patterns.
    status: pending
  - id: r-usage-doc
    content: Add a short usage note + minimal R snippet in the script header (and optionally README) for generating the boxplot from TOTAL rows.
    status: pending
---

## Goal (detailed execution spec)

Create a **single, Step1D-style, HPC-safe utility script** that generates an **R-ready “SNPs per sample” count table** from a folder of per-chromosome VCFs.

### What the script must do (no ambiguity)

- **Input**: a directory containing multiple compressed VCFs (`*.vcf.gz`) where the filename contains the substring `chr` (case-insensitive), e.g. `chr01.vcf.gz`, `Chr01_consolidated.vcf.gz`, `mycohort_chr3.filtered.vcf.gz`.
- **Detect**: all matching `*.vcf.gz` in that directory (non-recursive), sort them deterministically, and treat each file as a “chromosome-like unit” (we do not assume chromosomes 0–17; we follow filenames).
- **Count definition (fixed)**: for each detected VCF file, compute per-sample counts of **called genotypes** at sites that satisfy:
- biallelic SNP sites only (`-m2 -M2 -v snps`)
- `FILTER=PASS` only (`-f PASS`)
- The main “SNP count per sample” metric is:
- `n_called = nRefHom + nNonRefHom + nHets`
- i.e. number of non-missing genotypes among the filtered biallelic PASS SNP sites.
- **Output**: write **one single file** (not 17 files) that contains:
- Per-file (“per chr”) rows for each `(sample, chr_label)`
- Plus a per-sample **TOTAL** row summing across all `chr_label`s
- **R-friendliness**: the default output must be a **single gzip-compressed CSV** (`.csv.gz`) that can be read with `data.table::fread()` without additional parsing.
- **HPC safety**: it must run correctly when executed via `sbatch` from SLURM spool locations; it must not rely on `$PWD` being the repo root.

### Script name + location (fixed)

- **Create**: `GATK_Pipeline_KH_v1/utils/snp_counts_per_sample_from_chr_vcfs.sh`
- Must be runnable both ways:
- interactive: `bash snp_counts_per_sample_from_chr_vcfs.sh --vcf-dir /path/to/vcfs`
- SLURM: `sbatch --export=ALL,VCF_DIR=/path/to/vcfs ... snp_counts_per_sample_from_chr_vcfs.sh`

## Design constraints pulled from the pipeline’s HPC playbook

From [`docs/past_problems_and_resolutions.md`](/Users/khangha/Library/CloudStorage/OneDrive-TheUniversityofQueensland/Code/Apple_GATK_2025/KH/GATK_Pipeline_KH_v1/docs/past_problems_and_resolutions.md):

- Always support **absolute-path safety** for SLURM spool execution; rely on `${PIPELINE_ROOT}` when sourcing `lib/`/`config/`.
- Prefer **`module purge` then `module load bcftools/...`** to avoid “illegal instruction” or picking up wrong binaries.
- Avoid adding fragile Python dependencies for core QC; bcftools + awk is the safest baseline.

## “Step1d-style” conventions to mirror

From [`modules/step1d/templates/master_vcf_analysis.sh`](/Users/khangha/Library/CloudStorage/OneDrive-TheUniversityofQueensland/Code/Apple_GATK_2025/KH/GATK_Pipeline_KH_v1/modules/step1d/templates/master_vcf_analysis.sh):

- `set -Eeuo pipefail`, optional `--dry-run` mode, env-var-overridable config.
- Auto-detection patterns (it already falls back to `find` when pattern yields nothing).
- `module purge` then load `BCFTOOLS_MODULE` with default `bcftools/1.18-gcc-12.3.0`.
- Optional structured logging if `lib/logging.sh` is available.

## What the new script will do

### CLI + environment contract (exact)

Support both CLI flags and env vars; **CLI flags override env vars**; env vars override defaults.

- `--vcf-dir <dir>` (env: `VCF_DIR`, default: `$PWD`)
- Directory scanned for inputs (non-recursive).
- `--work-dir <dir>` (env: `WORK_DIR`, default: `VCF_DIR`)
- Where the output file is written.
- `--out <path>` (env: `OUT_FILE`, default: `${WORK_DIR}/snps_per_sample_counts.csv.gz`)
- Output path; must end in `.csv.gz` by default.
- `--name-filter <pattern>` (env: `VCF_NAME_FILTER`, default: `*chr*`)
- Case-insensitive substring-style filter for filenames (the default means “contains chr”).
- This is applied in addition to the `*.vcf.gz` suffix requirement.
- `--bcftools-module <module>` (env: `BCFTOOLS_MODULE`, default: `bcftools/1.18-gcc-12.3.0`)
- Exact string to pass to `module load`.
- `--samples <file>` (env: `SAMPLES_FILE`, default: empty)
- Optional: if provided, pass to `bcftools stats -s <file>` to restrict samples to a list.
- If not provided, use `-s -` (all samples).
- `--dry-run` (env: `DRY_RUN`, default: false)
- Print intended actions and exit 0 without creating output.
- `--help`

### Input discovery (exact)

Inside `VCF_DIR`, collect candidates using logic equivalent to:

- Candidates must satisfy:
- filename matches `*.vcf.gz`
- and filename matches the case-insensitive filter `*chr*` (default)
- Deterministic ordering:
- sort by filename with `LC_ALL=C sort`
- For each file:
- `source_vcf = absolute path`
- `chr_label = basename(source_vcf)` with `.vcf.gz` stripped (only that suffix; do not strip extra tokens)

If **no files** match, exit non-zero with a clear error.

### Inputs

- **Folder** containing per-chromosome (or chr-named) `*.vcf.gz` files.
- Defaults configurable via env vars (Step1d-style):
- `VCF_DIR` (default: `$PWD`)
- `VCF_NAME_FILTER` (default: case-insensitive match containing `chr`)
- `BCFTOOLS_MODULE` (default: `bcftools/1.18-gcc-12.3.0`)
- `OUT_FILE` (default: `${WORK_DIR}/snps_per_sample_counts.csv.gz`)
- `WORK_DIR` (default: `VCF_DIR`)
- `DRY_RUN` (default: false)

### Detection

- Discover inputs with something equivalent to:
- `find "$VCF_DIR" -maxdepth 1 -type f -name "*.vcf.gz" -iname "*chr*" | sort`
- Extract a `chr_label` for each file from its basename (strip `.vcf.gz`).

### Counting definition (your chosen default)

For each detected VCF file, the script must run a counting pipeline equivalent to:

- Filter stream:
- `bcftools view -m2 -M2 -v snps -f PASS -Ou "${vcf}"`
- Stats computation:
- pipe into `bcftools stats`:
- `bcftools stats -s <samples> -`
- Where `<samples>` is `-` (all samples) or `SAMPLES_FILE` (if provided)

**Important parsing requirement (robust across bcftools versions):**

- Do **not** hardcode PSC column numbers.
- Instead:
- Parse the `bcftools stats` output to find the header line that describes PSC fields (it’s typically a comment line starting with `# PSC` that enumerates columns like `[4]nRefHom`, `[5]nNonRefHom`, etc.).
- Build an index map for at least these PSC fields:
- `sample`
- `nRefHom`
- `nNonRefHom`
- `nHets`
- `nMissing` (or the equivalent missing count field if named slightly differently)
- Then parse `PSC` records using the discovered indices.

Compute metrics per `(sample, chr_label)` row:

- `n_refhom = nRefHom`
- `n_homalt = nNonRefHom`
- `n_het = nHets`
- `n_missing = nMissing`
- `n_called = n_refhom + n_homalt + n_het`
- `n_total = n_called + n_missing`

### Single output file (R-friendly)

Write **one** compressed CSV (easy for `data.table::fread`):

- **Per-chromosome rows**: one row per `(sample, chr_label)`
- **TOTAL rows**: one row per sample with `chr_label=TOTAL` summing across all chromosomes

Recommended columns:

- `scope` (`CHR` or `TOTAL`)
- `chr_label`
- `sample`
- `n_called` (the main boxplot variable)
- `n_missing`
- `n_total`
- `n_refhom`, `n_homalt`, `n_het` (useful QC breakdown)
- `source_vcf` (only for `scope=CHR`, blank for TOTAL)

**CSV details (exact):**

- Include a header row.
- Use comma delimiter.
- Quote only when necessary (standard CSV).
- Gzip-compress the final output so the final path is `*.csv.gz`.
- Deterministic row order:
- First all `scope=CHR` rows sorted by `chr_label`, then by `sample`
- Then all `scope=TOTAL` rows sorted by `sample`

## How it plugs into your R boxplot

- In R: `dt <- fread("snps_per_sample_counts.csv.gz"); dt_total <- dt[chr_label=="TOTAL"]; ggplot(dt_total, aes(y=n_called, x="")) + geom_boxplot() + geom_jitter(width=0.15, alpha=0.3)`

## Validation and guard rails

- Validate at runtime:
- bcftools exists after module load (`command -v bcftools`).
- detected file count > 0.
- sample sets are consistent across VCFs:
- get sample list from first file via `bcftools query -l`
- for each subsequent file, compare lists; if mismatch, exit non-zero and print which file differs
- indexing expectations:
- do not assume `.tbi` exists; bcftools can stream without random access
- but if an input file is corrupt/unreadable, fail fast on that file with a clear message
- Log the exact filters applied and number of files processed.

### Aggregation to TOTAL rows (exact)

After producing all CHR rows in memory/temporary file:

- For each `sample`, sum across all CHR rows:
- `n_called_total = sum(n_called)`
- `n_missing_total = sum(n_missing)`
- `n_total_total = sum(n_total)`
- and similarly sum `n_refhom`, `n_homalt`, `n_het`
- Emit one TOTAL row per sample with:
- `scope=TOTAL`
- `chr_label=TOTAL`
- `source_vcf=` (empty)

## Location

- Add the script under [`GATK_Pipeline_KH_v1/utils/`](/Users/khangha/Library/CloudStorage/OneDrive-TheUniversityofQueensland/Code/Apple_GATK_2025/KH/GATK_Pipeline_KH_v1/utils/) so it matches your existing utilities (`vcfstats.sh`, `bcfstats.sh`) and is easy to call from Step1d or wrappers later.

