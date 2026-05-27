#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --job-name=plink2_ld_decay
#SBATCH --time=1-00:00:00
#SBATCH --partition=general
#SBATCH --account=a_qaafi_cas
#SBATCH -o slurm.plink2_ld_decay.%j.out
#SBATCH -e slurm.plink2_ld_decay.%j.err

set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  # Recommended: one PLINK2 run against the concatenated Chr01-Chr17 VCF.
  sbatch utils/poplddecay.sh --combined-vcf <Chr01-Chr17.vcf.gz> --outdir <dir> [options]

  # Optional: process individual per-chromosome VCFs sequentially in one job.
  sbatch utils/poplddecay.sh --vcf-dir <dir> --outdir <dir> [options]

  # Optional: process one chromosome only.
  bash utils/poplddecay.sh --vcf <Chr01.vcf.gz> --outdir <dir> --prefix Chr01 [options]
  bash utils/poplddecay.sh --vcf-dir <dir> --chr Chr01 --outdir <dir> [options]

Purpose:
  Despite the historical filename, this script now computes LD decay with PLINK2,
  not PopLDdecay. The primary path consumes one concatenated Chr01-Chr17 VCF.
  Per-chromosome input paths remain available for focused or legacy runs.

Input selection:
  --combined-vcf PATH    One VCF/BCF containing Chr01-Chr17; processed once.
  --vcf PATH             Explicit one-chromosome VCF/BCF for a single run.
  --vcf-dir DIR          Directory containing per-chromosome VCFs.
  --vcf-pattern PATTERN  Filename pattern used with --vcf-dir.
                         Default: Chr%02d_merged_noCommon_snps.vcf.gz
                         The %02d token is replaced by chromosome number.
  --chr CHR              Run only one chromosome, e.g. Chr01 or 1.
  --chrom-start INT      First chromosome number for multi-chrom run. Default: 1
  --chrom-end INT        Last chromosome number for multi-chrom run. Default: 17

Output and reporting:
  --outdir DIR           Output directory. Default: ./plink2_ld_decay
  --prefix NAME          Output prefix. Only valid with one chromosome.
  --combined-prefix NAME Output prefix for --combined-vcf or aggregated
                         per-chromosome results. Default: plink2_ld_decay_combined
  --no-rds               Skip RDS output.
  --keep-pairwise        Copy large PLINK2 pairwise LD .vcor files to --outdir.
  --keep-tmp             Keep the per-run working directory for debugging.

PLINK2 LD parameters:
  --ld-window-kb INT     Maximum LD distance in kb. Default: 300
  --bin-kb INT           LD decay bin width in kb. Default: 10
  --maf FLOAT            PLINK2 --maf threshold. Default: 0.05
  --geno FLOAT           PLINK2 --geno variant missingness threshold. Default: 1
  --mind FLOAT           PLINK2 --mind sample missingness threshold. Default: 1
  --ld-window-r2 FLOAT   Minimum r2 included in pairwise report. Default: 0
  --thin FLOAT           Randomly keep this fraction of variants. Default: 0.1
  --threads INT          PLINK2 threads. Default: SLURM_CPUS_PER_TASK or 8

Environment:
  --conda-env NAME       Conda env containing modern plink2 and Rscript. Default: plink2_ld
  --miniforge-module MOD Module used to activate Conda. Default: miniforge/26.1.0-0
  --no-conda             Do not activate Conda; use plink2/Rscript already on PATH.

Outputs:
  <outdir>/<combined-prefix>.plink2_ld_decay.tsv  Table from --combined-vcf,
                                                   retaining chromosome labels
  <outdir>/<combined-prefix>.plink2_ld_decay.rds  RDS from --combined-vcf
  <outdir>/Chr01.plink2_ld_decay.tsv           Per-chrom binned LD table
  <outdir>/Chr01.plink2_ld_decay.rds           Per-chrom RDS list
  <outdir>/Chr01.plink2_ld_decay.log           Per-chrom command log
  <outdir>/Chr01.plink2_ld_decay.steps.tsv     Per-chrom audit report
  <outdir>/<combined-prefix>.tsv               Combined Chr01-Chr17 binned table
  <outdir>/<combined-prefix>.rds               Combined RDS list for QMD/R
  <outdir>/<combined-prefix>.steps.tsv         Whole-run audit report
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
  local report_file="$1"
  local chromosome="$2"
  local step="$3"
  local status="$4"
  local started_at="${5:-0}"
  local detail="${6:-}"
  local ended_at
  local duration=""
  local safe_detail

  ended_at="$(date +%s)"
  if [[ "${started_at}" =~ ^[0-9]+$ && "${started_at}" -gt 0 ]]; then
    duration="$((ended_at - started_at))"
  fi

  safe_detail="$(printf '%s' "${detail}" | tr '\t\n' '  ')"
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$(timestamp)" \
    "${chromosome}" \
    "${step}" \
    "${status}" \
    "${duration}" \
    "${safe_detail}" >> "${report_file}"
}

ensure_conda_tools() {
  if [[ "${USE_CONDA}" == "true" ]]; then
    if command -v module >/dev/null 2>&1; then
      module purge || true
      module load "${MINIFORGE_MODULE}" || true
    fi
    if [[ -n "${ROOTMINIFORGE:-}" && -f "${ROOTMINIFORGE}/etc/profile.d/conda.sh" ]]; then
      source "${ROOTMINIFORGE}/etc/profile.d/conda.sh"
      conda activate "${CONDA_ENV}"
    elif command -v conda >/dev/null 2>&1; then
      eval "$(conda shell.bash hook)"
      conda activate "${CONDA_ENV}"
    else
      echo "[plink2_ld_decay] ERROR: conda not found. Load Miniforge or use --no-conda." >&2
      exit 1
    fi
  fi

  if ! command -v plink2 >/dev/null 2>&1; then
    echo "[plink2_ld_decay] ERROR: plink2 not found on PATH." >&2
    exit 1
  fi
  if [[ "${WRITE_RDS}" == "true" ]] && ! command -v Rscript >/dev/null 2>&1; then
    echo "[plink2_ld_decay] ERROR: Rscript not found on PATH. Use --no-rds or activate an R environment." >&2
    exit 1
  fi
  if ! command -v rsync >/dev/null 2>&1; then
    echo "[plink2_ld_decay] ERROR: rsync not found on PATH." >&2
    exit 1
  fi
}

COMBINED_VCF=""
VCF=""
VCF_DIR=""
VCF_PATTERN='Chr%02d_merged_noCommon_snps.vcf.gz'
CHR_ARG=""
CHROM_START=1
CHROM_END=17
OUTDIR="plink2_ld_decay"
PREFIX=""
COMBINED_PREFIX="plink2_ld_decay_combined"
LD_WINDOW_KB=300
BIN_KB=10
MAF=0.05
GENO=1
MIND=1
LD_WINDOW_R2=0
THIN=0.1
THREADS="${SLURM_CPUS_PER_TASK:-8}"
CONDA_ENV="${PLINK2_LD_CONDA_ENV:-plink2_ld}"
MINIFORGE_MODULE="${MINIFORGE_MODULE:-miniforge/26.1.0-0}"
USE_CONDA=true
WRITE_RDS=true
KEEP_PAIRWISE=false
KEEP_TMP=false

while [[ $# -gt 0 ]]; do
  case "$1" in
    --combined-vcf) COMBINED_VCF="$2"; shift 2 ;;
    --vcf) VCF="$2"; shift 2 ;;
    --vcf-dir) VCF_DIR="$2"; shift 2 ;;
    --vcf-pattern) VCF_PATTERN="$2"; shift 2 ;;
    --chr) CHR_ARG="$2"; shift 2 ;;
    --chrom-start) CHROM_START="$2"; shift 2 ;;
    --chrom-end) CHROM_END="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --prefix) PREFIX="$2"; shift 2 ;;
    --combined-prefix) COMBINED_PREFIX="$2"; shift 2 ;;
    --ld-window-kb) LD_WINDOW_KB="$2"; shift 2 ;;
    --bin-kb) BIN_KB="$2"; shift 2 ;;
    --maf) MAF="$2"; shift 2 ;;
    --geno) GENO="$2"; shift 2 ;;
    --mind) MIND="$2"; shift 2 ;;
    --ld-window-r2) LD_WINDOW_R2="$2"; shift 2 ;;
    --thin) THIN="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --conda-env) CONDA_ENV="$2"; shift 2 ;;
    --miniforge-module) MINIFORGE_MODULE="$2"; shift 2 ;;
    --no-conda) USE_CONDA=false; shift ;;
    --no-rds) WRITE_RDS=false; shift ;;
    --keep-pairwise) KEEP_PAIRWISE=true; shift ;;
    --keep-tmp) KEEP_TMP=true; shift ;;
    --help|-h) usage; exit 0 ;;
    *) echo "[plink2_ld_decay] ERROR: unknown argument: $1" >&2; usage >&2; exit 1 ;;
  esac
done

if [[ -n "${COMBINED_VCF}" && ( -n "${VCF}" || -n "${VCF_DIR}" || -n "${CHR_ARG}" ) ]]; then
  echo "[plink2_ld_decay] ERROR: --combined-vcf cannot be combined with --vcf, --vcf-dir, or --chr." >&2
  exit 1
fi

if [[ -n "${PREFIX}" && -z "${VCF}" && -z "${CHR_ARG}" ]]; then
  echo "[plink2_ld_decay] ERROR: --prefix is only valid with a one-chromosome --vcf or --chr run; use --combined-prefix for --combined-vcf." >&2
  exit 1
fi

mkdir -p "${OUTDIR}"
FINAL_OUTDIR="$(cd "${OUTDIR}" && pwd)"
RUN_REPORT="${FINAL_OUTDIR%/}/${COMBINED_PREFIX}.steps.tsv"
printf 'timestamp\tchromosome\tstep\tstatus\tduration_seconds\tdetail\n' > "${RUN_REPORT}"

WORK_PARENT="${TMPDIR:-${FINAL_OUTDIR}}"
mkdir -p "${WORK_PARENT}"
WORK_DIR="$(mktemp -d "${WORK_PARENT%/}/plink2_ld_decay.${SLURM_JOB_ID:-$$}.XXXXXX")"

cleanup() {
  if [[ "${KEEP_TMP}" != "true" ]]; then
    rm -rf "${WORK_DIR}"
  else
    echo "[plink2_ld_decay] Keeping temporary directory: ${WORK_DIR}" >&2
  fi
}
trap cleanup EXIT

ensure_conda_tools
report_step "${RUN_REPORT}" "ALL" "initialise" "DONE" 0 "work_dir=${WORK_DIR}; outdir=${FINAL_OUTDIR}; plink2=$(command -v plink2); thin=${THIN}; ld_window_kb=${LD_WINDOW_KB}; bin_kb=${BIN_KB}"

CHR_NUMS=()
VCF_PATHS=()
CHR_LABELS=()
PREFIXES=()
INPUT_MODES=()

if [[ -n "${COMBINED_VCF}" ]]; then
  if [[ ! -f "${COMBINED_VCF}" ]]; then
    echo "[plink2_ld_decay] ERROR: combined input VCF not found: ${COMBINED_VCF}" >&2
    exit 1
  fi
  CHR_NUMS+=("0")
  VCF_PATHS+=("${COMBINED_VCF}")
  CHR_LABELS+=("ALL")
  PREFIXES+=("${COMBINED_PREFIX}")
  INPUT_MODES+=("combined")
elif [[ -n "${VCF}" ]]; then
  if [[ ! -f "${VCF}" ]]; then
    echo "[plink2_ld_decay] ERROR: input VCF not found: ${VCF}" >&2
    exit 1
  fi
  if [[ -n "${CHR_ARG}" ]]; then
    chr_num="$(normalise_chr_number "${CHR_ARG}")"
    chr_label="$(chr_label_from_number "${chr_num}")"
  else
    base_for_chr="$(basename "${VCF}")"
    chr_guess="$(printf '%s' "${base_for_chr}" | sed -nE 's/.*[Cc][Hh][Rr]0*([0-9]+).*/\1/p' | head -n 1)"
    if [[ -n "${chr_guess}" ]]; then
      chr_num="${chr_guess}"
      chr_label="$(chr_label_from_number "${chr_num}")"
    else
      chr_num="0"
      chr_label="ChrNA"
    fi
  fi
  CHR_NUMS+=("${chr_num}")
  VCF_PATHS+=("${VCF}")
  CHR_LABELS+=("${chr_label}")
  PREFIXES+=("${PREFIX:-${chr_label}}")
  INPUT_MODES+=("chromosome")
elif [[ -n "${CHR_ARG}" ]]; then
  if [[ -z "${VCF_DIR}" ]]; then
    echo "[plink2_ld_decay] ERROR: --vcf-dir is required with --chr." >&2
    exit 1
  fi
  chr_num="$(normalise_chr_number "${CHR_ARG}")"
  chr_label="$(chr_label_from_number "${chr_num}")"
  vcf_basename="$(printf "${VCF_PATTERN}" "${chr_num}")"
  vcf_path="${VCF_DIR%/}/${vcf_basename}"
  CHR_NUMS+=("${chr_num}")
  VCF_PATHS+=("${vcf_path}")
  CHR_LABELS+=("${chr_label}")
  PREFIXES+=("${PREFIX:-${chr_label}}")
  INPUT_MODES+=("chromosome")
else
  if [[ -z "${VCF_DIR}" ]]; then
    echo "[plink2_ld_decay] ERROR: provide --combined-vcf, --vcf, or --vcf-dir." >&2
    usage >&2
    exit 1
  fi
  for chr_num in $(seq "${CHROM_START}" "${CHROM_END}"); do
    chr_label="$(chr_label_from_number "${chr_num}")"
    vcf_basename="$(printf "${VCF_PATTERN}" "${chr_num}")"
    vcf_path="${VCF_DIR%/}/${vcf_basename}"
    CHR_NUMS+=("${chr_num}")
    VCF_PATHS+=("${vcf_path}")
    CHR_LABELS+=("${chr_label}")
    PREFIXES+=("${chr_label}")
    INPUT_MODES+=("chromosome")
  done
fi

SUMMARY_FILES=()
RDS_FILES=()

for i in "${!VCF_PATHS[@]}"; do
  chr_num="${CHR_NUMS[$i]}"
  chr_label="${CHR_LABELS[$i]}"
  vcf_path="${VCF_PATHS[$i]}"
  sample_prefix="${PREFIXES[$i]}"
  input_mode="${INPUT_MODES[$i]}"

  if [[ ! -f "${vcf_path}" ]]; then
    echo "[plink2_ld_decay] ERROR: input VCF not found for ${chr_label}: ${vcf_path}" >&2
    report_step "${RUN_REPORT}" "${chr_label}" "check_input" "FAILED" 0 "missing=${vcf_path}"
    exit 1
  fi

  chr_report="${FINAL_OUTDIR%/}/${sample_prefix}.plink2_ld_decay.steps.tsv"
  log_file="${FINAL_OUTDIR%/}/${sample_prefix}.plink2_ld_decay.log"
  printf 'timestamp\tchromosome\tstep\tstatus\tduration_seconds\tdetail\n' > "${chr_report}"

  work_prefix="${WORK_DIR%/}/${sample_prefix}.plink2_ld"
  summary_unsorted="${WORK_DIR%/}/${sample_prefix}.plink2_ld_decay.unsorted.tsv"
  summary_tsv="${WORK_DIR%/}/${sample_prefix}.plink2_ld_decay.tsv"
  rds_file="${WORK_DIR%/}/${sample_prefix}.plink2_ld_decay.rds"

  report_step "${chr_report}" "${chr_label}" "initialise" "DONE" 0 "input=${vcf_path}; work_dir=${WORK_DIR}; outdir=${FINAL_OUTDIR}"
  report_step "${RUN_REPORT}" "${chr_label}" "chromosome_start" "START" 0 "input=${vcf_path}"

  {
    echo "[plink2_ld_decay] Input: ${vcf_path}"
    echo "[plink2_ld_decay] Input mode: ${input_mode}"
    echo "[plink2_ld_decay] Chromosome scope: ${chr_label}"
    echo "[plink2_ld_decay] Work directory: ${WORK_DIR}"
    echo "[plink2_ld_decay] Final output prefix: ${FINAL_OUTDIR%/}/${sample_prefix}"
    echo "[plink2_ld_decay] PLINK2: $(command -v plink2)"
  } | tee "${log_file}"

  step_start="$(date +%s)"
  report_step "${chr_report}" "${chr_label}" "run_plink2_ld" "START" 0 "output_prefix=${work_prefix}"
  plink2 \
    --vcf "${vcf_path}" \
    --double-id \
    --allow-extra-chr \
    --set-all-var-ids '@:#:$r:$a' \
    --snps-only just-acgt \
    --max-alleles 2 \
    --maf "${MAF}" \
    --geno "${GENO}" \
    --mind "${MIND}" \
    --thin "${THIN}" \
    --r2-unphased \
    --ld-window-kb "${LD_WINDOW_KB}" \
    --ld-window-r2 "${LD_WINDOW_R2}" \
    --threads "${THREADS}" \
    --out "${work_prefix}" 2>&1 | tee -a "${log_file}"
  report_step "${chr_report}" "${chr_label}" "run_plink2_ld" "DONE" "${step_start}" "output_prefix=${work_prefix}"

  pairwise_ld=""
  for candidate in "${work_prefix}.vcor" "${work_prefix}"*.vcor "${work_prefix}.vcor.zst" "${work_prefix}"*.vcor.zst; do
    if [[ -f "${candidate}" ]]; then
      pairwise_ld="${candidate}"
      break
    fi
  done
  if [[ -z "${pairwise_ld}" ]]; then
    echo "[plink2_ld_decay] ERROR: Could not find PLINK2 .vcor output for ${work_prefix}" | tee -a "${log_file}"
    report_step "${chr_report}" "${chr_label}" "find_pairwise_ld" "FAILED" 0 "prefix=${work_prefix}"
    exit 1
  fi

  step_start="$(date +%s)"
  report_step "${chr_report}" "${chr_label}" "summarise_ld_decay" "START" 0 "pairwise_ld=${pairwise_ld}; bin_kb=${BIN_KB}"
  if [[ "${pairwise_ld}" == *.zst ]]; then
    reader=(zstdcat "${pairwise_ld}")
  else
    reader=(cat "${pairwise_ld}")
  fi

  "${reader[@]}" | awk -v bin_bp="$((BIN_KB * 1000))" -v fixed_chr="${chr_label}" -v input_mode="${input_mode}" '
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
      chrom_a = col["CHROM_A"];
      chrom_b = col["CHROM_B"];
      if (input_mode == "combined" && (chrom_a == "" || chrom_b == "")) {
        print "Missing CHROM_A or CHROM_B columns in multi-chromosome PLINK2 output" > "/dev/stderr";
        exit 2;
      }
      next;
    }
    {
      r2 = $r2_col;
      if (r2 == "nan" || r2 == "NA" || r2 == ".") next;
      out_chr = fixed_chr;
      if (input_mode == "combined") {
        if ($chrom_a != $chrom_b) next;
        out_chr = $chrom_a;
      }
      dist = $pos_b - $pos_a;
      if (dist < 0) dist = -dist;
      bin = int(dist / bin_bp) * bin_bp;
      key = out_chr SUBSEP bin;
      sum[key] += r2;
      n[key] += 1;
    }
    END {
      for (key in n) {
        split(key, parts, SUBSEP);
        bin = parts[2];
        print parts[1], bin, bin + bin_bp - 1, (bin + (bin_bp / 2)) / 1000, n[key], sum[key] / n[key];
      }
    }
  ' > "${summary_unsorted}"

  {
    head -n 1 "${summary_unsorted}"
    tail -n +2 "${summary_unsorted}" | sort -k1,1 -k2,2n
  } > "${summary_tsv}"
  report_step "${chr_report}" "${chr_label}" "summarise_ld_decay" "DONE" "${step_start}" "summary=${summary_tsv}"

  if [[ "${WRITE_RDS}" == "true" ]]; then
    step_start="$(date +%s)"
    report_step "${chr_report}" "${chr_label}" "write_rds" "START" 0 "summary=${summary_tsv}; rds=${rds_file}"
    Rscript - "${summary_tsv}" "${rds_file}" "${chr_label}" "${sample_prefix}" "${vcf_path}" "${LD_WINDOW_KB}" "${BIN_KB}" "${MAF}" "${GENO}" "${MIND}" "${LD_WINDOW_R2}" "${THIN}" <<'RSCRIPT'
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
thin <- as.numeric(args[[12]])
ld_decay <- read.delim(summary_tsv, check.names = FALSE)
metadata <- list(
  created_at = as.character(Sys.time()), chromosome = chromosome,
  sample_prefix = sample_prefix, input_vcf = input_vcf,
  ld_window_kb = ld_window_kb, bin_kb = bin_kb, maf = maf, geno = geno,
  mind = mind, ld_window_r2 = ld_window_r2, thin = thin,
  source_tsv = basename(summary_tsv), n_rows = nrow(ld_decay)
)
saveRDS(list(ld_decay = ld_decay, metadata = metadata), rds_file)
RSCRIPT
    report_step "${chr_report}" "${chr_label}" "write_rds" "DONE" "${step_start}" "rds=${rds_file}"
  else
    report_step "${chr_report}" "${chr_label}" "write_rds" "SKIP" 0 "--no-rds set"
  fi

  step_start="$(date +%s)"
  report_step "${chr_report}" "${chr_label}" "copy_outputs" "START" 0 "destination=${FINAL_OUTDIR}; keep_pairwise=${KEEP_PAIRWISE}"
  rsync -rhivPt "${summary_tsv}" "${FINAL_OUTDIR}/" 2>&1 | tee -a "${log_file}"
  if [[ "${WRITE_RDS}" == "true" ]]; then
    rsync -rhivPt "${rds_file}" "${FINAL_OUTDIR}/" 2>&1 | tee -a "${log_file}"
    RDS_FILES+=("${FINAL_OUTDIR%/}/$(basename "${rds_file}")")
  fi
  if [[ "${KEEP_PAIRWISE}" == "true" ]]; then
    rsync -rhivPt "${pairwise_ld}" "${FINAL_OUTDIR}/" 2>&1 | tee -a "${log_file}"
  fi
  report_step "${chr_report}" "${chr_label}" "copy_outputs" "DONE" "${step_start}" "destination=${FINAL_OUTDIR}"
  report_step "${chr_report}" "${chr_label}" "complete" "DONE" 0 "summary=${FINAL_OUTDIR%/}/$(basename "${summary_tsv}")"
  report_step "${RUN_REPORT}" "${chr_label}" "chromosome_complete" "DONE" 0 "summary=${FINAL_OUTDIR%/}/$(basename "${summary_tsv}")"
  SUMMARY_FILES+=("${FINAL_OUTDIR%/}/$(basename "${summary_tsv}")")

  echo "[plink2_ld_decay] Completed ${chr_label}: ${FINAL_OUTDIR%/}/$(basename "${summary_tsv}")" | tee -a "${log_file}"
done

if [[ "${#SUMMARY_FILES[@]}" -gt 1 ]]; then
  combined_tsv="${WORK_DIR%/}/${COMBINED_PREFIX}.tsv"
  {
    head -n 1 "${SUMMARY_FILES[0]}"
    for f in "${SUMMARY_FILES[@]}"; do
      tail -n +2 "${f}"
    done
  } > "${combined_tsv}"
  rsync -rhivPt "${combined_tsv}" "${FINAL_OUTDIR}/" >/dev/null

  if [[ "${WRITE_RDS}" == "true" ]]; then
    combined_rds="${WORK_DIR%/}/${COMBINED_PREFIX}.rds"
    Rscript - "${combined_tsv}" "${combined_rds}" "${VCF_DIR}" "${LD_WINDOW_KB}" "${BIN_KB}" "${MAF}" "${GENO}" "${MIND}" "${LD_WINDOW_R2}" "${THIN}" <<'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)
combined_tsv <- args[[1]]
combined_rds <- args[[2]]
vcf_dir <- args[[3]]
ld_window_kb <- as.numeric(args[[4]])
bin_kb <- as.numeric(args[[5]])
maf <- as.numeric(args[[6]])
geno <- as.numeric(args[[7]])
mind <- as.numeric(args[[8]])
ld_window_r2 <- as.numeric(args[[9]])
thin <- as.numeric(args[[10]])
ld_decay <- read.delim(combined_tsv, check.names = FALSE)
ld_decay <- ld_decay[order(ld_decay$chromosome, ld_decay$bin_start_bp), , drop = FALSE]
metadata <- list(
  created_at = as.character(Sys.time()), source_dir = vcf_dir,
  chromosomes = sort(unique(ld_decay$chromosome)), n_rows = nrow(ld_decay),
  ld_window_kb = ld_window_kb, bin_kb = bin_kb, maf = maf, geno = geno,
  mind = mind, ld_window_r2 = ld_window_r2, thin = thin,
  source_tsv = basename(combined_tsv)
)
saveRDS(list(ld_decay = ld_decay, metadata = metadata), combined_rds)
RSCRIPT
    rsync -rhivPt "${combined_rds}" "${FINAL_OUTDIR}/" >/dev/null
  fi
fi

report_step "${RUN_REPORT}" "ALL" "complete" "DONE" 0 "outdir=${FINAL_OUTDIR}; chromosomes=${#SUMMARY_FILES[@]}"
echo "[plink2_ld_decay] All done. Results: ${FINAL_OUTDIR}"
