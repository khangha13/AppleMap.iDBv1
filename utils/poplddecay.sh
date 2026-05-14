#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --job-name=plink2_ld_decay
#SBATCH --time=1-00:00:00
#SBATCH --partition=general
#SBATCH --account=a_qaafi_cas
#SBATCH --array=1-17
#SBATCH -o slurm.plink2_ld_decay.%A_%a.out
#SBATCH -e slurm.plink2_ld_decay.%A_%a.err

set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  # Recommended Slurm array use, one chromosome per task:
  sbatch utils/poplddecay.sh --vcf-dir <dir> --outdir <dir> [options]

  # Single chromosome / interactive use:
  bash utils/poplddecay.sh --vcf <Chr01.vcf.gz> --outdir <dir> --prefix Chr01 [options]
  bash utils/poplddecay.sh --vcf-dir <dir> --chr Chr01 --outdir <dir> [options]

Purpose:
  Despite the historical filename, this script now computes LD decay with PLINK2,
  not PopLDdecay. It is designed for conventional per-chromosome LD decay runs.
  Each chromosome task produces a small binned TSV suitable for R/QMD plotting.

Input selection:
  --vcf PATH             Explicit input VCF/BCF for a single run.
  --vcf-dir DIR          Directory containing per-chromosome VCFs.
  --vcf-pattern PATTERN  Filename pattern used with --vcf-dir.
                         Default: Chr%02d_merged_noCommon_snps.vcf.gz
                         The %02d token is replaced by chromosome number.
  --chr CHR              Chromosome label for single run, e.g. Chr01 or 1.
                         In Slurm array mode, SLURM_ARRAY_TASK_ID is used.
  --chrom-start INT      First chromosome number. Default: 1
  --chrom-end INT        Last chromosome number. Default: 17

Output and reporting:
  --outdir DIR           Output directory. Default: ./plink2_ld_decay
  --prefix NAME          Output prefix. Default: chromosome label or VCF basename.
  --no-rds               Skip per-chromosome RDS output.
  --combine-rds          Combine existing per-chromosome TSVs in --outdir into
                         a single all-chromosome RDS, then exit.
  --combined-rds NAME    Filename for --combine-rds. Default:
                         plink2_ld_decay_combined.rds
  --keep-pairwise        Copy the large PLINK2 pairwise LD .vcor file to --outdir.
  --keep-tmp             Keep the per-task working directory for debugging.

PLINK2 LD parameters:
  --ld-window-kb INT     Maximum LD distance in kb. Default: 100
  --bin-kb INT           LD decay bin width in kb. Default: 10
  --maf FLOAT            PLINK2 --maf threshold. Default: 0.05
  --geno FLOAT           PLINK2 --geno variant missingness threshold. Default: 1
  --mind FLOAT           PLINK2 --mind sample missingness threshold. Default: 1
  --ld-window-r2 FLOAT   Minimum r2 included in pairwise report. Default: 0
  --threads INT          PLINK2 threads. Default: SLURM_CPUS_PER_TASK or 8

Environment:
  --plink-module MOD     PLINK2 module. Default: plink/2.00a3.6-gcc-11.3.0
  --miniforge-module MOD Module used to activate R env. Default: miniforge/26.1.0-0
  --r-env NAME           Conda env containing Rscript. Default: rplot
  --no-r-env             Do not activate Conda before writing RDS; use Rscript on PATH.

Outputs per chromosome:
  <outdir>/<prefix>.plink2_ld_decay.tsv        Binned LD decay table
  <outdir>/<prefix>.plink2_ld_decay.rds        Binned LD decay RDS list
  <outdir>/<prefix>.plink2_ld_decay.log        Command log
  <outdir>/<prefix>.plink2_ld_decay.steps.tsv  Per-step audit report
  <outdir>/<prefix>.*.vcor.gz                  Optional with --keep-pairwise

QMD plotting:
  Run --combine-rds after the array finishes, then read the combined RDS in QMD.
EOF
}

timestamp() {
  date '+%Y-%m-%dT%H:%M:%S%z'
}

normalise_chr_number() {
  local chr="$1"
  chr="${chr#Chr}"
  chr="${chr#chr}"
  printf '%d' "${chr}"
}

chr_label_from_number() {
  local chr_num="$1"
  printf 'Chr%02d' "${chr_num}"
}

report_step() {
  local step="$1"
  local status="$2"
  local started_at="${3:-0}"
  local detail="${4:-}"
  local ended_at
  local duration=""
  local safe_detail

  ended_at="$(date +%s)"
  if [[ "${started_at}" =~ ^[0-9]+$ && "${started_at}" -gt 0 ]]; then
    duration="$((ended_at - started_at))"
  fi

  safe_detail="$(printf '%s' "${detail}" | tr '\t\n' '  ')"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$(timestamp)" \
    "${CHR_LABEL:-NA}" \
    "${sample_prefix:-global}" \
    "${step}" \
    "${status}" \
    "${duration}" \
    "${safe_detail}" >> "${STEP_REPORT}"
}

ensure_rscript() {
  if command -v Rscript >/dev/null 2>&1; then
    return 0
  fi

  if [[ "${USE_R_ENV}" == "true" ]]; then
    if command -v module >/dev/null 2>&1; then
      module load "${MINIFORGE_MODULE}" || true
    fi
    if [[ -n "${ROOTMINIFORGE:-}" && -f "${ROOTMINIFORGE}/etc/profile.d/conda.sh" ]]; then
      source "${ROOTMINIFORGE}/etc/profile.d/conda.sh"
      conda activate "${R_ENV}"
    elif command -v conda >/dev/null 2>&1; then
      eval "$(conda shell.bash hook)"
      conda activate "${R_ENV}"
    fi
  fi

  if ! command -v Rscript >/dev/null 2>&1; then
    echo "[plink2_ld_decay] ERROR: Rscript not found. Use --no-rds, or load/activate an R environment." >&2
    exit 1
  fi
}

VCF=""
VCF_DIR=""
VCF_PATTERN='Chr%02d_merged_noCommon_snps.vcf.gz'
CHR_ARG=""
CHROM_START=1
CHROM_END=17
OUTDIR="plink2_ld_decay"
PREFIX=""
LD_WINDOW_KB=100
BIN_KB=10
MAF=0.05
GENO=1
MIND=1
LD_WINDOW_R2=0
THREADS="${SLURM_CPUS_PER_TASK:-8}"
PLINK_MODULE="${PLINK_MODULE:-plink/2.00a3.6-gcc-11.3.0}"
MINIFORGE_MODULE="${MINIFORGE_MODULE:-miniforge/26.1.0-0}"
R_ENV="${R_ENV:-rplot}"
USE_R_ENV=true
WRITE_RDS=true
COMBINE_RDS=false
COMBINED_RDS_NAME="plink2_ld_decay_combined.rds"
KEEP_PAIRWISE=false
KEEP_TMP=false

while [[ $# -gt 0 ]]; do
  case "$1" in
    --vcf)
      VCF="$2"
      shift 2
      ;;
    --vcf-dir)
      VCF_DIR="$2"
      shift 2
      ;;
    --vcf-pattern)
      VCF_PATTERN="$2"
      shift 2
      ;;
    --chr)
      CHR_ARG="$2"
      shift 2
      ;;
    --chrom-start)
      CHROM_START="$2"
      shift 2
      ;;
    --chrom-end)
      CHROM_END="$2"
      shift 2
      ;;
    --outdir)
      OUTDIR="$2"
      shift 2
      ;;
    --prefix)
      PREFIX="$2"
      shift 2
      ;;
    --ld-window-kb)
      LD_WINDOW_KB="$2"
      shift 2
      ;;
    --bin-kb)
      BIN_KB="$2"
      shift 2
      ;;
    --maf)
      MAF="$2"
      shift 2
      ;;
    --geno)
      GENO="$2"
      shift 2
      ;;
    --mind)
      MIND="$2"
      shift 2
      ;;
    --ld-window-r2)
      LD_WINDOW_R2="$2"
      shift 2
      ;;
    --threads)
      THREADS="$2"
      shift 2
      ;;
    --plink-module)
      PLINK_MODULE="$2"
      shift 2
      ;;
    --miniforge-module)
      MINIFORGE_MODULE="$2"
      shift 2
      ;;
    --r-env)
      R_ENV="$2"
      shift 2
      ;;
    --no-r-env)
      USE_R_ENV=false
      shift
      ;;
    --no-rds)
      WRITE_RDS=false
      shift
      ;;
    --combine-rds)
      COMBINE_RDS=true
      shift
      ;;
    --combined-rds)
      COMBINED_RDS_NAME="$2"
      shift 2
      ;;
    --keep-pairwise)
      KEEP_PAIRWISE=true
      shift
      ;;
    --keep-tmp)
      KEEP_TMP=true
      shift
      ;;
    --help|-h)
      usage
      exit 0
      ;;
    *)
      echo "[plink2_ld_decay] ERROR: unknown argument: $1" >&2
      usage >&2
      exit 1
      ;;
  esac
done

if [[ "${COMBINE_RDS}" == "true" ]]; then
  mkdir -p "${OUTDIR}"
  FINAL_OUTDIR="$(cd "${OUTDIR}" && pwd)"
  ensure_rscript
  Rscript - "${FINAL_OUTDIR}" "${COMBINED_RDS_NAME}" <<'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)
outdir <- args[[1]]
combined_name <- args[[2]]

files <- list.files(
  outdir,
  pattern = "^Chr[0-9]+\\.plink2_ld_decay\\.tsv$",
  full.names = TRUE
)
if (length(files) == 0) {
  stop("No per-chromosome TSV files found in ", outdir)
}

read_one <- function(path) {
  dat <- read.delim(path, check.names = FALSE)
  dat$source_file <- basename(path)
  dat
}

ld_decay <- do.call(rbind, lapply(files, read_one))
ld_decay <- ld_decay[order(ld_decay$chromosome, ld_decay$bin_start_bp), , drop = FALSE]

metadata <- list(
  created_at = as.character(Sys.time()),
  source_dir = outdir,
  source_files = basename(files),
  n_rows = nrow(ld_decay),
  chromosomes = sort(unique(ld_decay$chromosome))
)

saveRDS(
  list(ld_decay = ld_decay, metadata = metadata),
  file = file.path(outdir, combined_name)
)
message("Wrote ", file.path(outdir, combined_name))
RSCRIPT
  exit 0
fi

if [[ -z "${VCF}" ]]; then
  if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    CHR_NUM="${SLURM_ARRAY_TASK_ID}"
  elif [[ -n "${CHR_ARG}" ]]; then
    CHR_NUM="$(normalise_chr_number "${CHR_ARG}")"
  else
    echo "[plink2_ld_decay] ERROR: provide --vcf, --chr, or submit as a Slurm array task." >&2
    usage >&2
    exit 1
  fi

  if (( CHR_NUM < CHROM_START || CHR_NUM > CHROM_END )); then
    echo "[plink2_ld_decay] ERROR: chromosome ${CHR_NUM} outside ${CHROM_START}-${CHROM_END}." >&2
    exit 1
  fi

  if [[ -z "${VCF_DIR}" ]]; then
    echo "[plink2_ld_decay] ERROR: --vcf-dir is required unless --vcf is provided." >&2
    exit 1
  fi

  CHR_LABEL="$(chr_label_from_number "${CHR_NUM}")"
  vcf_basename="$(printf "${VCF_PATTERN}" "${CHR_NUM}")"
  VCF="${VCF_DIR%/}/${vcf_basename}"
else
  if [[ -n "${CHR_ARG}" ]]; then
    CHR_NUM="$(normalise_chr_number "${CHR_ARG}")"
    CHR_LABEL="$(chr_label_from_number "${CHR_NUM}")"
  else
    vcf_base_for_chr="$(basename "${VCF}")"
    chr_guess="$(printf '%s' "${vcf_base_for_chr}" | sed -nE 's/.*[Cc][Hh][Rr]0*([0-9]+).*/\1/p' | head -n 1)"
    if [[ -n "${chr_guess}" ]]; then
      CHR_NUM="${chr_guess}"
      CHR_LABEL="$(chr_label_from_number "${CHR_NUM}")"
    else
      CHR_NUM="NA"
      CHR_LABEL="NA"
    fi
  fi
fi

if [[ ! -f "${VCF}" ]]; then
  echo "[plink2_ld_decay] ERROR: input VCF not found: ${VCF}" >&2
  exit 1
fi

if [[ -z "${PREFIX}" ]]; then
  if [[ "${CHR_LABEL}" != "NA" ]]; then
    sample_prefix="${CHR_LABEL}"
  else
    base="$(basename "${VCF}")"
    sample_prefix="${base%.bcf}"
    sample_prefix="${sample_prefix%.vcf.gz}"
    sample_prefix="${sample_prefix%.vcf}"
  fi
else
  sample_prefix="${PREFIX}"
fi

if command -v module >/dev/null 2>&1; then
  module purge || true
  module load "${PLINK_MODULE}" || true
fi

if ! command -v plink2 >/dev/null 2>&1; then
  echo "[plink2_ld_decay] ERROR: plink2 not found on PATH. Load ${PLINK_MODULE} or set --plink-module." >&2
  exit 1
fi

if ! command -v rsync >/dev/null 2>&1; then
  echo "[plink2_ld_decay] ERROR: rsync not found on PATH." >&2
  exit 1
fi

mkdir -p "${OUTDIR}"
FINAL_OUTDIR="$(cd "${OUTDIR}" && pwd)"

WORK_PARENT="${TMPDIR:-${FINAL_OUTDIR}}"
mkdir -p "${WORK_PARENT}"
WORK_DIR="$(mktemp -d "${WORK_PARENT%/}/plink2_ld_decay.${SLURM_JOB_ID:-$$}.${SLURM_ARRAY_TASK_ID:-0}.XXXXXX")"

cleanup() {
  if [[ "${KEEP_TMP}" != "true" ]]; then
    rm -rf "${WORK_DIR}"
  else
    echo "[plink2_ld_decay] Keeping temporary directory: ${WORK_DIR}" >&2
  fi
}
trap cleanup EXIT

log_file="${FINAL_OUTDIR%/}/${sample_prefix}.plink2_ld_decay.log"
STEP_REPORT="${FINAL_OUTDIR%/}/${sample_prefix}.plink2_ld_decay.steps.tsv"
summary_unsorted="${WORK_DIR%/}/${sample_prefix}.plink2_ld_decay.unsorted.tsv"
summary_tsv="${WORK_DIR%/}/${sample_prefix}.plink2_ld_decay.tsv"
rds_file="${WORK_DIR%/}/${sample_prefix}.plink2_ld_decay.rds"
work_prefix="${WORK_DIR%/}/${sample_prefix}.plink2_ld"

printf 'timestamp\tchromosome\tsample_prefix\tstep\tstatus\tduration_seconds\tdetail\n' > "${STEP_REPORT}"
report_step "initialise" "DONE" 0 "input=${VCF}; work_dir=${WORK_DIR}; outdir=${FINAL_OUTDIR}; ld_window_kb=${LD_WINDOW_KB}; bin_kb=${BIN_KB}; maf=${MAF}; geno=${GENO}; mind=${MIND}; ld_window_r2=${LD_WINDOW_R2}; threads=${THREADS}"

{
  echo "[plink2_ld_decay] Input: ${VCF}"
  echo "[plink2_ld_decay] Chromosome: ${CHR_LABEL}"
  echo "[plink2_ld_decay] Work directory: ${WORK_DIR}"
  echo "[plink2_ld_decay] Final output prefix: ${FINAL_OUTDIR%/}/${sample_prefix}"
  echo "[plink2_ld_decay] PLINK2: $(command -v plink2)"
} | tee "${log_file}"

step_start="$(date +%s)"
report_step "run_plink2_ld" "START" 0 "output_prefix=${work_prefix}"
plink2 \
  --vcf "${VCF}" \
  --double-id \
  --allow-extra-chr \
  --set-all-var-ids '@:#:$r:$a' \
  --snps-only just-acgt \
  --max-alleles 2 \
  --maf "${MAF}" \
  --geno "${GENO}" \
  --mind "${MIND}" \
  --r2-unphased gz \
  --ld-window-kb "${LD_WINDOW_KB}" \
  --ld-window-r2 "${LD_WINDOW_R2}" \
  --bad-ld \
  --threads "${THREADS}" \
  --out "${work_prefix}" 2>&1 | tee -a "${log_file}"
report_step "run_plink2_ld" "DONE" "${step_start}" "output_prefix=${work_prefix}"

pairwise_ld=""
for candidate in \
  "${work_prefix}.vcor.gz" \
  "${work_prefix}.unphased.vcor.gz" \
  "${work_prefix}"*.vcor.gz \
  "${work_prefix}.vcor" \
  "${work_prefix}.unphased.vcor" \
  "${work_prefix}"*.vcor; do
  if [[ -f "${candidate}" ]]; then
    pairwise_ld="${candidate}"
    break
  fi
done

if [[ -z "${pairwise_ld}" ]]; then
  echo "[plink2_ld_decay] ERROR: Could not find PLINK2 .vcor output for ${work_prefix}" | tee -a "${log_file}"
  exit 1
fi

step_start="$(date +%s)"
report_step "summarise_ld_decay" "START" 0 "pairwise_ld=${pairwise_ld}; bin_kb=${BIN_KB}"
if [[ "${pairwise_ld}" == *.gz ]]; then
  reader=(gzip -cd "${pairwise_ld}")
else
  reader=(cat "${pairwise_ld}")
fi

"${reader[@]}" | awk -v bin_bp="$((BIN_KB * 1000))" -v chr="${CHR_LABEL}" '
  BEGIN {
    OFS = "\t";
    print "chromosome", "bin_start_bp", "bin_end_bp", "bin_mid_kb", "n_pairs", "mean_r2";
  }
  NR == 1 {
    for (i = 1; i <= NF; i++) {
      name = $i;
      sub(/^#/, "", name);
      col[name] = i;
    }
    pos_a = col["POS_A"];
    pos_b = col["POS_B"];
    r2_col = col["UNPHASED_R2"];
    if (pos_a == "" || pos_b == "" || r2_col == "") {
      print "Missing POS_A, POS_B, or UNPHASED_R2 columns in PLINK2 output" > "/dev/stderr";
      exit 2;
    }
    next;
  }
  {
    r2 = $r2_col;
    if (r2 == "nan" || r2 == "NA" || r2 == ".") next;
    dist = $pos_b - $pos_a;
    if (dist < 0) dist = -dist;
    bin = int(dist / bin_bp) * bin_bp;
    sum[bin] += r2;
    n[bin] += 1;
  }
  END {
    for (bin in n) {
      print chr, bin, bin + bin_bp - 1, (bin + (bin_bp / 2)) / 1000, n[bin], sum[bin] / n[bin];
    }
  }
' > "${summary_unsorted}"

{
  head -n 1 "${summary_unsorted}"
  tail -n +2 "${summary_unsorted}" | sort -k2,2n
} > "${summary_tsv}"
report_step "summarise_ld_decay" "DONE" "${step_start}" "summary=${summary_tsv}"

if [[ "${WRITE_RDS}" == "true" ]]; then
  step_start="$(date +%s)"
  report_step "write_rds" "START" 0 "summary=${summary_tsv}; rds=${rds_file}"
  ensure_rscript
  Rscript - "${summary_tsv}" "${rds_file}" "${CHR_LABEL}" "${sample_prefix}" "${VCF}" "${LD_WINDOW_KB}" "${BIN_KB}" "${MAF}" "${GENO}" "${MIND}" "${LD_WINDOW_R2}" <<'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)
summary_tsv <- args[[1]]
rds_file <- args[[2]]
chromosome <- args[[3]]
sample_prefix <- args[[4]]
input_vcf <- args[[5]]
ld_window_kb <- as.numeric(args[[6]])
bin_kb <- as.numeric(args[[7]])
maf <- as.numeric(args[[8]])
geno <- as.numeric(args[[9]])
mind <- as.numeric(args[[10]])
ld_window_r2 <- as.numeric(args[[11]])

ld_decay <- read.delim(summary_tsv, check.names = FALSE)
metadata <- list(
  created_at = as.character(Sys.time()),
  chromosome = chromosome,
  sample_prefix = sample_prefix,
  input_vcf = input_vcf,
  ld_window_kb = ld_window_kb,
  bin_kb = bin_kb,
  maf = maf,
  geno = geno,
  mind = mind,
  ld_window_r2 = ld_window_r2,
  source_tsv = basename(summary_tsv),
  n_rows = nrow(ld_decay)
)

saveRDS(list(ld_decay = ld_decay, metadata = metadata), rds_file)
RSCRIPT
  report_step "write_rds" "DONE" "${step_start}" "rds=${rds_file}"
else
  report_step "write_rds" "SKIP" 0 "--no-rds set"
fi

step_start="$(date +%s)"
report_step "copy_outputs" "START" 0 "destination=${FINAL_OUTDIR}; keep_pairwise=${KEEP_PAIRWISE}"
echo "[plink2_ld_decay] Copying outputs back to ${FINAL_OUTDIR}" | tee -a "${log_file}"
rsync -rhivPt "${summary_tsv}" "${FINAL_OUTDIR}/" 2>&1 | tee -a "${log_file}"
if [[ "${WRITE_RDS}" == "true" ]]; then
  rsync -rhivPt "${rds_file}" "${FINAL_OUTDIR}/" 2>&1 | tee -a "${log_file}"
fi
if [[ "${KEEP_PAIRWISE}" == "true" ]]; then
  rsync -rhivPt "${pairwise_ld}" "${FINAL_OUTDIR}/" 2>&1 | tee -a "${log_file}"
fi
report_step "copy_outputs" "DONE" "${step_start}" "destination=${FINAL_OUTDIR}"

report_step "complete" "DONE" 0 "summary=${FINAL_OUTDIR%/}/${sample_prefix}.plink2_ld_decay.tsv; rds=${FINAL_OUTDIR%/}/${sample_prefix}.plink2_ld_decay.rds; log=${log_file}; step_report=${STEP_REPORT}"
echo "[plink2_ld_decay] Done: ${FINAL_OUTDIR%/}/${sample_prefix}.plink2_ld_decay.tsv" | tee -a "${log_file}"
