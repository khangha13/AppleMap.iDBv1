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
  --combined-prefix NAME Output prefix for --combined-vcf.
                         Default: plink2_ld_decay_combined
  --keep-tmp             Keep the per-run working directory for debugging.

PLINK2 LD parameters:
  --ld-window-kb INT     Maximum LD distance in kb. Default: 300
  --maf FLOAT            PLINK2 --maf threshold. Default: 0.05
  --geno FLOAT           PLINK2 --geno variant missingness threshold. Default: 1
  --mind FLOAT           PLINK2 --mind sample missingness threshold. Default: 1
  --ld-window-r2 FLOAT   Minimum r2 included in pairwise report. Default: 0
  --thin FLOAT           Randomly keep this fraction of variants. Default: 0.1
  --threads INT          PLINK2 threads. Default: SLURM_CPUS_PER_TASK or 8

Environment:
  --conda-env NAME       Conda env containing modern plink2. Default: plink2_ld
  --miniforge-module MOD Module used to activate Conda. Default: miniforge/26.1.0-0
  --no-conda             Do not activate Conda; use plink2 already on PATH.

Outputs:
  <outdir>/<combined-prefix>.plink2_ld.vcor[.zst]  Raw LD from --combined-vcf
  <outdir>/Chr01.plink2_ld.vcor[.zst]              Raw LD from a chromosome run

  A multi-chromosome --vcf-dir run writes one raw .vcor file per chromosome.
  Use --combined-vcf when one graphing input file is required.
  Under sbatch, progress and PLINK2 console output are captured in the Slurm
  .out file configured at the top of this script.
EOF
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
MAF=0.05
GENO=1
MIND=1
LD_WINDOW_R2=0
THIN=0.1
THREADS="${SLURM_CPUS_PER_TASK:-8}"
CONDA_ENV="${PLINK2_LD_CONDA_ENV:-plink2_ld}"
MINIFORGE_MODULE="${MINIFORGE_MODULE:-miniforge/26.1.0-0}"
USE_CONDA=true
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
    --maf) MAF="$2"; shift 2 ;;
    --geno) GENO="$2"; shift 2 ;;
    --mind) MIND="$2"; shift 2 ;;
    --ld-window-r2) LD_WINDOW_R2="$2"; shift 2 ;;
    --thin) THIN="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --conda-env) CONDA_ENV="$2"; shift 2 ;;
    --miniforge-module) MINIFORGE_MODULE="$2"; shift 2 ;;
    --no-conda) USE_CONDA=false; shift ;;
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

VCF_PATHS=()
CHR_LABELS=()
PREFIXES=()

if [[ -n "${COMBINED_VCF}" ]]; then
  if [[ ! -f "${COMBINED_VCF}" ]]; then
    echo "[plink2_ld_decay] ERROR: combined input VCF not found: ${COMBINED_VCF}" >&2
    exit 1
  fi
  VCF_PATHS+=("${COMBINED_VCF}")
  CHR_LABELS+=("ALL")
  PREFIXES+=("${COMBINED_PREFIX}")
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
  VCF_PATHS+=("${VCF}")
  CHR_LABELS+=("${chr_label}")
  PREFIXES+=("${PREFIX:-${chr_label}}")
elif [[ -n "${CHR_ARG}" ]]; then
  if [[ -z "${VCF_DIR}" ]]; then
    echo "[plink2_ld_decay] ERROR: --vcf-dir is required with --chr." >&2
    exit 1
  fi
  chr_num="$(normalise_chr_number "${CHR_ARG}")"
  chr_label="$(chr_label_from_number "${chr_num}")"
  vcf_basename="$(printf "${VCF_PATTERN}" "${chr_num}")"
  vcf_path="${VCF_DIR%/}/${vcf_basename}"
  VCF_PATHS+=("${vcf_path}")
  CHR_LABELS+=("${chr_label}")
  PREFIXES+=("${PREFIX:-${chr_label}}")
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
    VCF_PATHS+=("${vcf_path}")
    CHR_LABELS+=("${chr_label}")
    PREFIXES+=("${chr_label}")
  done
fi

for i in "${!VCF_PATHS[@]}"; do
  chr_label="${CHR_LABELS[$i]}"
  vcf_path="${VCF_PATHS[$i]}"
  sample_prefix="${PREFIXES[$i]}"

  if [[ ! -f "${vcf_path}" ]]; then
    echo "[plink2_ld_decay] ERROR: input VCF not found for ${chr_label}: ${vcf_path}" >&2
    exit 1
  fi

  work_prefix="${WORK_DIR%/}/${sample_prefix}.plink2_ld"

  echo "[plink2_ld_decay] Input: ${vcf_path}"
  echo "[plink2_ld_decay] Chromosome scope: ${chr_label}"
  echo "[plink2_ld_decay] Work directory: ${WORK_DIR}"
  echo "[plink2_ld_decay] Final output prefix: ${FINAL_OUTDIR%/}/${sample_prefix}.plink2_ld"
  echo "[plink2_ld_decay] PLINK2: $(command -v plink2)"

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
    --out "${work_prefix}"

  pairwise_ld=""
  for candidate in "${work_prefix}.vcor" "${work_prefix}"*.vcor "${work_prefix}.vcor.zst" "${work_prefix}"*.vcor.zst; do
    if [[ -f "${candidate}" ]]; then
      pairwise_ld="${candidate}"
      break
    fi
  done
  if [[ -z "${pairwise_ld}" ]]; then
    echo "[plink2_ld_decay] ERROR: Could not find PLINK2 .vcor output for ${work_prefix}" >&2
    exit 1
  fi

  final_pairwise="${FINAL_OUTDIR%/}/$(basename "${pairwise_ld}")"
  mv -f "${pairwise_ld}" "${final_pairwise}"
  echo "[plink2_ld_decay] Completed ${chr_label}: ${final_pairwise}"
done

echo "[plink2_ld_decay] All done. Results: ${FINAL_OUTDIR}"
