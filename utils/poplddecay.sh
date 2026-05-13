#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --job-name=poplddecay
#SBATCH --time=1-00:00:00
#SBATCH --partition=general
#SBATCH --account=a_qaafi_chs
#SBATCH -o slurm.poplddecay.%j.out
#SBATCH -e slurm.poplddecay.%j.err

set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  bash utils/poplddecay.sh --vcf <input.vcf.gz> [options]
  bash utils/poplddecay.sh <input1.vcf.gz> [input2.vcf.gz ...]

Options:
  --vcf PATH             Input VCF/BCF. Can be provided multiple times.
  --outdir DIR           Output directory. Default: ./ld_decay
  --prefix NAME          Output prefix. Only valid with one input VCF.
  --max-dist INT         PopLDdecay -MaxDist value. Default: 1000
  --maf FLOAT            Optional PopLDdecay -MAF value.
  --miss FLOAT           Optional PopLDdecay -Miss value.
  --no-filter            Do not pre-filter to biallelic SNPs with bcftools.
  --conda-env NAME       Conda env containing PopLDdecay. Default: poplddecay
  --no-conda             Do not activate Conda; use commands already on PATH.
  --miniforge-module MOD Module to load before conda activate. Default:
                         miniforge/26.1.0-0
  --bcftools-module MOD  Optional bcftools module. Default:
                         bcftools/1.18-gcc-12.3.0
  --help                 Show this help.

Outputs per input:
  <outdir>/<prefix>.filtered_snps.vcf.gz       temporary filtered VCF
  <outdir>/<prefix>.poplddecay.stat.gz         PopLDdecay statistics
  <outdir>/<prefix>.ld_decay.*                 Plot_OnePop.pl outputs
  <outdir>/<prefix>.poplddecay.log             Command log

Notes:
  - Default filtering keeps biallelic SNPs only:
      bcftools view -m2 -M2 -v snps
  - If bcftools is unavailable, use --no-filter and provide a SNP-only VCF.
EOF
}

OUTDIR="ld_decay"
PREFIX=""
MAX_DIST="1000"
MAF=""
MISS=""
FILTER_SNPS=true
USE_CONDA=true
CONDA_ENV="${POPLDDECAY_CONDA_ENV:-poplddecay}"
MINIFORGE_MODULE="${MINIFORGE_MODULE:-miniforge/26.1.0-0}"
BCFTOOLS_MODULE="${BCFTOOLS_MODULE:-bcftools/1.18-gcc-12.3.0}"
VCFS=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --vcf)
      VCFS+=("$2")
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
    --max-dist)
      MAX_DIST="$2"
      shift 2
      ;;
    --maf)
      MAF="$2"
      shift 2
      ;;
    --miss)
      MISS="$2"
      shift 2
      ;;
    --no-filter)
      FILTER_SNPS=false
      shift
      ;;
    --conda-env)
      CONDA_ENV="$2"
      shift 2
      ;;
    --no-conda)
      USE_CONDA=false
      shift
      ;;
    --miniforge-module)
      MINIFORGE_MODULE="$2"
      shift 2
      ;;
    --bcftools-module)
      BCFTOOLS_MODULE="$2"
      shift 2
      ;;
    --help|-h)
      usage
      exit 0
      ;;
    -*)
      echo "[poplddecay] ERROR: unknown option: $1" >&2
      usage >&2
      exit 1
      ;;
    *)
      VCFS+=("$1")
      shift
      ;;
  esac
done

if [[ "${#VCFS[@]}" -eq 0 ]]; then
  echo "[poplddecay] ERROR: provide at least one VCF." >&2
  usage >&2
  exit 1
fi

if [[ -n "${PREFIX}" && "${#VCFS[@]}" -gt 1 ]]; then
  echo "[poplddecay] ERROR: --prefix is only valid with one input VCF." >&2
  exit 1
fi

for vcf in "${VCFS[@]}"; do
  if [[ ! -f "${vcf}" ]]; then
    echo "[poplddecay] ERROR: input VCF not found: ${vcf}" >&2
    exit 1
  fi
done

if command -v module >/dev/null 2>&1; then
  module purge || true
  if [[ "${USE_CONDA}" == "true" ]]; then
    module load "${MINIFORGE_MODULE}" || true
  fi
  if [[ "${FILTER_SNPS}" == "true" ]]; then
    module load "${BCFTOOLS_MODULE}" || true
  fi
fi

if [[ "${USE_CONDA}" == "true" ]]; then
  if [[ -n "${ROOTMINIFORGE:-}" && -f "${ROOTMINIFORGE}/etc/profile.d/conda.sh" ]]; then
    source "${ROOTMINIFORGE}/etc/profile.d/conda.sh"
    conda activate "${CONDA_ENV}"
  elif command -v conda >/dev/null 2>&1; then
    eval "$(conda shell.bash hook)"
    conda activate "${CONDA_ENV}"
  else
    echo "[poplddecay] ERROR: conda not found. Load Miniforge, or rerun with --no-conda if PopLDdecay is already on PATH." >&2
    exit 1
  fi
fi

if ! command -v PopLDdecay >/dev/null 2>&1; then
  echo "[poplddecay] ERROR: PopLDdecay not found on PATH." >&2
  exit 1
fi

if ! command -v Plot_OnePop.pl >/dev/null 2>&1; then
  echo "[poplddecay] ERROR: Plot_OnePop.pl not found on PATH." >&2
  exit 1
fi

if [[ "${FILTER_SNPS}" == "true" ]] && ! command -v bcftools >/dev/null 2>&1; then
  echo "[poplddecay] ERROR: bcftools not found. Load bcftools or use --no-filter." >&2
  exit 1
fi

mkdir -p "${OUTDIR}"

for vcf in "${VCFS[@]}"; do
  if [[ -n "${PREFIX}" ]]; then
    sample_prefix="${PREFIX}"
  else
    base="$(basename "${vcf}")"
    sample_prefix="${base%.bcf}"
    sample_prefix="${sample_prefix%.vcf.gz}"
    sample_prefix="${sample_prefix%.vcf}"
  fi

  log_file="${OUTDIR%/}/${sample_prefix}.poplddecay.log"
  stat_file="${OUTDIR%/}/${sample_prefix}.poplddecay.stat.gz"
  plot_prefix="${OUTDIR%/}/${sample_prefix}.ld_decay"

  echo "[poplddecay] Input: ${vcf}" | tee "${log_file}"
  echo "[poplddecay] Output prefix: ${OUTDIR%/}/${sample_prefix}" | tee -a "${log_file}"

  run_vcf="${vcf}"
  if [[ "${FILTER_SNPS}" == "true" ]]; then
    filtered_vcf="${OUTDIR%/}/${sample_prefix}.filtered_snps.vcf.gz"
    echo "[poplddecay] Filtering to biallelic SNPs: ${filtered_vcf}" | tee -a "${log_file}"
    bcftools view -m2 -M2 -v snps -Oz -o "${filtered_vcf}" "${vcf}" 2>&1 | tee -a "${log_file}"
    bcftools index --tbi --force "${filtered_vcf}" 2>&1 | tee -a "${log_file}"
    run_vcf="${filtered_vcf}"
  fi

  cmd=(PopLDdecay -InVCF "${run_vcf}" -OutStat "${stat_file}" -MaxDist "${MAX_DIST}")
  if [[ -n "${MAF}" ]]; then
    cmd+=(-MAF "${MAF}")
  fi
  if [[ -n "${MISS}" ]]; then
    cmd+=(-Miss "${MISS}")
  fi

  echo "[poplddecay] Running: ${cmd[*]}" | tee -a "${log_file}"
  "${cmd[@]}" 2>&1 | tee -a "${log_file}"

  echo "[poplddecay] Plotting: Plot_OnePop.pl -inFile ${stat_file} -output ${plot_prefix}" | tee -a "${log_file}"
  Plot_OnePop.pl -inFile "${stat_file}" -output "${plot_prefix}" 2>&1 | tee -a "${log_file}"

  echo "[poplddecay] Done: ${stat_file}" | tee -a "${log_file}"
done
