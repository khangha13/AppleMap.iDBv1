#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=256G
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
  --keep-vcor            Retained for compatibility; raw PLINK2 LD is copied
                         to --outdir by default.
  --keep-tmp             Keep the per-run working directory for debugging.

PLINK2 LD parameters:
  --ld-window-kb INT     Full-distance PLINK2 maximum LD distance in kb.
                         Default: 1000
  --fine-ld-window-kb INT
                         Fine-profile PLINK2 maximum LD distance in kb.
                         Default: 5
  --bin-bp INT           Width of the full-distance profile bins in bp.
                         Default: 100. Fine profiles are always 50 bp,
                         10 bp, and 5 bp over 0-5 kb.
  --maf FLOAT            PLINK2 --maf threshold. Default: 0.05
  --geno FLOAT           PLINK2 --geno variant missingness threshold. Default: 1
  --mind FLOAT           PLINK2 --mind sample missingness threshold. Default: 1
  --ld-window-r2 FLOAT   Minimum r2 included in pairwise report. Default: 0
  --thin FLOAT           Full-distance run: randomly keep this fraction of
                         variants. Default: 0.2
  --fine-thin FLOAT      Fine-profile run: randomly keep this fraction of
                         variants. Default: 1
  --threads INT          PLINK2 threads. Default: SLURM_CPUS_PER_TASK or 8

Environment:
  --conda-env NAME       Conda env containing modern plink2. Default: plink2_ld
  --miniforge-module MOD Module used to activate Conda. Default: miniforge/26.1.0-0
  --no-conda             Do not activate Conda; use plink2 already on PATH.

Outputs:
  <outdir>/<prefix>.plink2_ld_decay.tsv           Long-format binned LD table.
                         Columns: profile, bin_bp, max_distance_bp, chromosome,
                         bin_start_bp, bin_end_bp, bin_mid_kb, n_pairs, mean_r2.
                         Profiles:
                           full_<bin-bp>bp_<ld-window-kb>kb
                           fine_50bp_<fine-ld-window-kb>kb
                           fine_10bp_<fine-ld-window-kb>kb
                           fine_5bp_<fine-ld-window-kb>kb
  <outdir>/<prefix>.plink2_ld.vcor.zst            Raw PLINK2 LD copied by default

  A multi-chromosome --vcf-dir run writes one long-format summary table per chromosome.
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
LD_WINDOW_KB=1000
FINE_LD_WINDOW_KB=5
BIN_BP=100
MAF=0.05
GENO=1
MIND=1
LD_WINDOW_R2=0
THIN=0.2
FINE_THIN=1
THREADS="${SLURM_CPUS_PER_TASK:-8}"
CONDA_ENV="${PLINK2_LD_CONDA_ENV:-plink2_ld}"
MINIFORGE_MODULE="${MINIFORGE_MODULE:-miniforge/26.1.0-0}"
USE_CONDA=true
KEEP_VCOR=true
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
    --fine-ld-window-kb) FINE_LD_WINDOW_KB="$2"; shift 2 ;;
    --bin-bp) BIN_BP="$2"; shift 2 ;;
    --maf) MAF="$2"; shift 2 ;;
    --geno) GENO="$2"; shift 2 ;;
    --mind) MIND="$2"; shift 2 ;;
    --ld-window-r2) LD_WINDOW_R2="$2"; shift 2 ;;
    --thin) THIN="$2"; shift 2 ;;
    --fine-thin) FINE_THIN="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --conda-env) CONDA_ENV="$2"; shift 2 ;;
    --miniforge-module) MINIFORGE_MODULE="$2"; shift 2 ;;
    --no-conda) USE_CONDA=false; shift ;;
    --keep-vcor) KEEP_VCOR=true; shift ;;
    --keep-tmp) KEEP_TMP=true; shift ;;
    --help|-h) usage; exit 0 ;;
    *) echo "[plink2_ld_decay] ERROR: unknown argument: $1" >&2; usage >&2; exit 1 ;;
  esac
done

if [[ ! "${BIN_BP}" =~ ^[1-9][0-9]*$ ]]; then
  echo "[plink2_ld_decay] ERROR: --bin-bp must be a positive integer number of base pairs." >&2
  exit 1
fi
if [[ ! "${LD_WINDOW_KB}" =~ ^[1-9][0-9]*$ ]]; then
  echo "[plink2_ld_decay] ERROR: --ld-window-kb must be a positive integer number of kilobases." >&2
  exit 1
fi
if [[ ! "${FINE_LD_WINDOW_KB}" =~ ^[1-9][0-9]*$ ]]; then
  echo "[plink2_ld_decay] ERROR: --fine-ld-window-kb must be a positive integer number of kilobases." >&2
  exit 1
fi

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

sort_ld_summary() {
  local unsorted="$1"
  local sorted="$2"
  local final="$3"

  {
    head -n 1 "${unsorted}"
    tail -n +2 "${unsorted}" | sort -k1,1 -k4,4 -k5,5n
  } > "${sorted}"
  mv -f "${sorted}" "${final}"
}

find_pairwise_ld() {
  local work_prefix="$1"
  local candidate

  for candidate in "${work_prefix}.vcor" "${work_prefix}"*.vcor "${work_prefix}.vcor.zst" "${work_prefix}"*.vcor.zst; do
    if [[ -f "${candidate}" ]]; then
      printf '%s\n' "${candidate}"
      return 0
    fi
  done
  return 1
}

append_ld_summary() {
  local pairwise_ld="$1"
  local fixed_chr="$2"
  local input_mode="$3"
  local summary_out="$4"
  local profile_mode="$5"
  local bin_bp="$6"
  local max_distance_bp="$7"
  local -a reader

  if [[ "${pairwise_ld}" == *.zst ]]; then
    reader=(plink2 --zst-decompress "${pairwise_ld}")
  else
    reader=(cat "${pairwise_ld}")
  fi

  "${reader[@]}" |
    awk \
      -v profile_mode="${profile_mode}" \
      -v bin_bp="${bin_bp}" \
      -v max_distance_bp="${max_distance_bp}" \
      -v fixed_chr="${fixed_chr}" \
      -v input_mode="${input_mode}" \
      -v main_out="${summary_out}" '
      BEGIN {
        OFS = "\t";
      }
      function record_bin(profile, out_chr, bin_width, profile_max_distance_bp, dist, r2) {
        bin = int(dist / bin_width) * bin_width;
        key = profile SUBSEP bin_width SUBSEP profile_max_distance_bp SUBSEP out_chr SUBSEP bin;
        sum[key] += r2;
        n[key] += 1;
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
        if (profile_mode == "full") {
          record_bin("full_" bin_bp "bp_" (max_distance_bp / 1000) "kb", out_chr, bin_bp, max_distance_bp, dist, r2);
        } else if (profile_mode == "fine" && dist < max_distance_bp) {
          record_bin("fine_50bp_" (max_distance_bp / 1000) "kb", out_chr, 50, max_distance_bp, dist, r2);
          record_bin("fine_10bp_" (max_distance_bp / 1000) "kb", out_chr, 10, max_distance_bp, dist, r2);
          record_bin("fine_5bp_" (max_distance_bp / 1000) "kb", out_chr, 5, max_distance_bp, dist, r2);
        }
      }
      END {
        for (key in n) {
          split(key, parts, SUBSEP);
          profile = parts[1];
          bin_width = parts[2] + 0;
          profile_max_distance_bp = parts[3] + 0;
          out_chr = parts[4];
          bin = parts[5] + 0;
          print profile, bin_width, profile_max_distance_bp, out_chr, bin, bin + bin_width - 1, (bin + bin_width / 2) / 1000, n[key], sum[key] / n[key] >> main_out;
        }
      }
    '
}

ensure_conda_tools

VCF_PATHS=()
CHR_LABELS=()
PREFIXES=()
INPUT_MODES=()

if [[ -n "${COMBINED_VCF}" ]]; then
  if [[ ! -f "${COMBINED_VCF}" ]]; then
    echo "[plink2_ld_decay] ERROR: combined input VCF not found: ${COMBINED_VCF}" >&2
    exit 1
  fi
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
    VCF_PATHS+=("${vcf_path}")
    CHR_LABELS+=("${chr_label}")
    PREFIXES+=("${chr_label}")
    INPUT_MODES+=("chromosome")
  done
fi

for i in "${!VCF_PATHS[@]}"; do
  chr_label="${CHR_LABELS[$i]}"
  vcf_path="${VCF_PATHS[$i]}"
  sample_prefix="${PREFIXES[$i]}"
  input_mode="${INPUT_MODES[$i]}"

  if [[ ! -f "${vcf_path}" ]]; then
    echo "[plink2_ld_decay] ERROR: input VCF not found for ${chr_label}: ${vcf_path}" >&2
    exit 1
  fi

  full_work_prefix="${WORK_DIR%/}/${sample_prefix}.plink2_ld.full"
  fine_work_prefix="${WORK_DIR%/}/${sample_prefix}.plink2_ld.fine"
  summary_unsorted="${WORK_DIR%/}/${sample_prefix}.plink2_ld_decay.unsorted.tsv"
  summary_tsv="${WORK_DIR%/}/${sample_prefix}.plink2_ld_decay.tsv"
  final_summary="${FINAL_OUTDIR%/}/${sample_prefix}.plink2_ld_decay.tsv"
  full_max_distance_bp=$((LD_WINDOW_KB * 1000))
  fine_max_distance_bp=$((FINE_LD_WINDOW_KB * 1000))

  echo "[plink2_ld_decay] Input: ${vcf_path}"
  echo "[plink2_ld_decay] Chromosome scope: ${chr_label}"
  echo "[plink2_ld_decay] Work directory: ${WORK_DIR}"
  echo "[plink2_ld_decay] Long-format summary output: ${final_summary}"
  echo "[plink2_ld_decay] PLINK2: $(command -v plink2)"

  printf 'profile\tbin_bp\tmax_distance_bp\tchromosome\tbin_start_bp\tbin_end_bp\tbin_mid_kb\tn_pairs\tmean_r2\n' > "${summary_unsorted}"
  full_thin_args=()
  fine_thin_args=()
  if [[ "${THIN}" != "1" && "${THIN}" != "1.0" ]]; then
    full_thin_args=(--thin "${THIN}")
  fi
  if [[ "${FINE_THIN}" != "1" && "${FINE_THIN}" != "1.0" ]]; then
    fine_thin_args=(--thin "${FINE_THIN}")
  fi

  echo "[plink2_ld_decay] Full-distance run: thin=${THIN}, window=${LD_WINDOW_KB}kb, bin=${BIN_BP}bp"
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
    "${full_thin_args[@]}" \
    --r2-unphased zs \
    --ld-window-kb "${LD_WINDOW_KB}" \
    --ld-window-r2 "${LD_WINDOW_R2}" \
    --threads "${THREADS}" \
    --out "${full_work_prefix}"

  if ! full_pairwise_ld="$(find_pairwise_ld "${full_work_prefix}")"; then
    echo "[plink2_ld_decay] ERROR: Could not find full-distance PLINK2 .vcor output for ${full_work_prefix}" >&2
    exit 1
  fi
  append_ld_summary "${full_pairwise_ld}" "${chr_label}" "${input_mode}" "${summary_unsorted}" "full" "${BIN_BP}" "${full_max_distance_bp}"

  echo "[plink2_ld_decay] Fine-profile run: thin=${FINE_THIN}, window=${FINE_LD_WINDOW_KB}kb, bins=50/10/5bp"
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
    "${fine_thin_args[@]}" \
    --r2-unphased zs \
    --ld-window-kb "${FINE_LD_WINDOW_KB}" \
    --ld-window-r2 "${LD_WINDOW_R2}" \
    --threads "${THREADS}" \
    --out "${fine_work_prefix}"

  if ! fine_pairwise_ld="$(find_pairwise_ld "${fine_work_prefix}")"; then
    echo "[plink2_ld_decay] ERROR: Could not find fine-profile PLINK2 .vcor output for ${fine_work_prefix}" >&2
    exit 1
  fi
  append_ld_summary "${fine_pairwise_ld}" "${chr_label}" "${input_mode}" "${summary_unsorted}" "fine" "${BIN_BP}" "${fine_max_distance_bp}"

  sort_ld_summary "${summary_unsorted}" "${summary_tsv}" "${final_summary}"

  if [[ "${KEEP_VCOR}" == "true" ]]; then
    for pairwise_ld in "${full_pairwise_ld}" "${fine_pairwise_ld}"; do
      final_pairwise="${FINAL_OUTDIR%/}/$(basename "${pairwise_ld}")"
      cp -f "${pairwise_ld}" "${final_pairwise}"
      echo "[plink2_ld_decay] Copied raw LD: ${final_pairwise}"
    done
  fi

  echo "[plink2_ld_decay] Completed ${chr_label}: ${final_summary}"
done

echo "[plink2_ld_decay] All done. Results: ${FINAL_OUTDIR}"
