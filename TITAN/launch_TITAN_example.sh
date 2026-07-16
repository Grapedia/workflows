#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
  cat <<'EOF'
Usage:
  ./launch_TITAN_example.sh [options] [-- extra_nextflow_options]

Required inputs can be provided either as options or environment variables:
  --output-dir PATH              TITAN_OUTPUT_DIR
  --previous-assembly PATH       TITAN_PREVIOUS_ASSEMBLY
  --new-assembly PATH            TITAN_NEW_ASSEMBLY
  --previous-annotations PATH    TITAN_PREVIOUS_ANNOTATIONS
  --rnaseq-samplesheet PATH      TITAN_RNASEQ_SAMPLESHEET
  --rnaseq-data-dir PATH         TITAN_RNASEQ_DATA_DIR
  --protein-samplesheet PATH     TITAN_PROTEIN_SAMPLESHEET
  --egapx-paramfile PATH         TITAN_EGAPX_PARAMFILE

Production-oriented options:
  --profile NAME                 Nextflow profile(s), for example slurm,apptainer
                                 Default: TITAN_PROFILE or slurm,apptainer
  --work-dir PATH                Nextflow work directory
                                 Default: TITAN_WORK_DIR or <output-dir>/work
  --run-name NAME                Stable Nextflow run name and report prefix
                                 Default: TITAN_RUN_NAME or titan_<timestamp>
  --module NAME                  Environment module to load before running Nextflow
                                 Default: TITAN_NEXTFLOW_MODULE or nextflow/24.04.3
                                 Set to empty string to disable module loading.
  --prepare-egapx-cache          Download EGAPx support data before launching TITAN
                                 Default: TITAN_PREPARE_EGAPX_CACHE or false
  --egapx-cache-dir PATH         Persistent EGAPx local cache directory
                                 Default: TITAN_EGAPX_CACHE_DIR or <project-dir>/.egapx_cache
  --egapx-runner-dir PATH        Local EGAPx runner source directory
                                 Default: TITAN_EGAPX_RUNNER_DIR or <project-dir>/.egapx_runner
  --resume                       Add -resume
  --force                        Allow writing into a non-empty output directory without --resume
  --dry-run                      Print the command instead of executing it
  -h, --help                     Show this help

Notes:
  TITAN has one public execution contract. Do not pass --workflow; evidence
  generation and Aegis are always launched together.

Example:
  TITAN_OUTPUT_DIR=/scratch/project/titan/riesling_hap2 \
  TITAN_PREVIOUS_ASSEMBLY=/data/ref/T2T_ref.fasta \
  TITAN_NEW_ASSEMBLY=/data/assemblies/riesling.hap2.fa \
  TITAN_PREVIOUS_ANNOTATIONS=/data/annotations/ref.gff3 \
  TITAN_RNASEQ_SAMPLESHEET=/data/rnaseq/RNAseq_samplesheet.csv \
  TITAN_RNASEQ_DATA_DIR=/data/rnaseq \
  TITAN_PROTEIN_SAMPLESHEET=/data/proteins/samplesheet.csv \
  TITAN_EGAPX_PARAMFILE=/data/input_egapx.yaml \
  ./launch_TITAN_example.sh --profile slurm,apptainer --resume
EOF
}

die() {
  echo "ERROR: $*" >&2
  exit 1
}

require_file() {
  local label="$1"
  local path="$2"
  [[ -n "$path" ]] || die "${label} is required"
  [[ -f "$path" ]] || die "${label} does not exist or is not a file: ${path}"
}

require_dir() {
  local label="$1"
  local path="$2"
  [[ -n "$path" ]] || die "${label} is required"
  [[ -d "$path" ]] || die "${label} does not exist or is not a directory: ${path}"
}

abs_path() {
  local path="$1"
  if [[ -d "$path" ]]; then
    cd "$path" && pwd -P
  else
    local dir
    dir="$(cd "$(dirname "$path")" && pwd -P)"
    printf '%s/%s\n' "$dir" "$(basename "$path")"
  fi
}

quote_cmd() {
  printf '%q ' "$@"
  printf '\n'
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
PROJECT_DIR="${TITAN_PROJECT_DIR:-$SCRIPT_DIR}"

PROFILE="${TITAN_PROFILE:-slurm,apptainer}"
OUTPUT_DIR="${TITAN_OUTPUT_DIR:-}"
PREVIOUS_ASSEMBLY="${TITAN_PREVIOUS_ASSEMBLY:-}"
NEW_ASSEMBLY="${TITAN_NEW_ASSEMBLY:-}"
PREVIOUS_ANNOTATIONS="${TITAN_PREVIOUS_ANNOTATIONS:-}"
RNASEQ_SAMPLESHEET="${TITAN_RNASEQ_SAMPLESHEET:-}"
RNASEQ_DATA_DIR="${TITAN_RNASEQ_DATA_DIR:-}"
PROTEIN_SAMPLESHEET="${TITAN_PROTEIN_SAMPLESHEET:-}"
EGAPX_PARAMFILE="${TITAN_EGAPX_PARAMFILE:-}"
RUN_NAME="${TITAN_RUN_NAME:-titan_$(date +%Y%m%d_%H%M%S)}"
NEXTFLOW_MODULE="${TITAN_NEXTFLOW_MODULE-nextflow/24.04.3}"
WORK_DIR="${TITAN_WORK_DIR:-}"
PREPARE_EGAPX_CACHE="${TITAN_PREPARE_EGAPX_CACHE:-false}"
EGAPX_CACHE_DIR="${TITAN_EGAPX_CACHE_DIR:-}"
EGAPX_RUNNER_DIR="${TITAN_EGAPX_RUNNER_DIR:-}"
RESUME=false
FORCE=false
DRY_RUN=false
EXTRA_NF_ARGS=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --output-dir) OUTPUT_DIR="${2:-}"; shift 2 ;;
    --previous-assembly) PREVIOUS_ASSEMBLY="${2:-}"; shift 2 ;;
    --new-assembly) NEW_ASSEMBLY="${2:-}"; shift 2 ;;
    --previous-annotations) PREVIOUS_ANNOTATIONS="${2:-}"; shift 2 ;;
    --rnaseq-samplesheet) RNASEQ_SAMPLESHEET="${2:-}"; shift 2 ;;
    --rnaseq-data-dir) RNASEQ_DATA_DIR="${2:-}"; shift 2 ;;
    --protein-samplesheet) PROTEIN_SAMPLESHEET="${2:-}"; shift 2 ;;
    --egapx-paramfile) EGAPX_PARAMFILE="${2:-}"; shift 2 ;;
    --profile) PROFILE="${2:-}"; shift 2 ;;
    --work-dir) WORK_DIR="${2:-}"; shift 2 ;;
    --run-name) RUN_NAME="${2:-}"; shift 2 ;;
    --module) NEXTFLOW_MODULE="${2:-}"; shift 2 ;;
    --prepare-egapx-cache) PREPARE_EGAPX_CACHE=true; shift ;;
    --egapx-cache-dir) EGAPX_CACHE_DIR="${2:-}"; shift 2 ;;
    --egapx-runner-dir) EGAPX_RUNNER_DIR="${2:-}"; shift 2 ;;
    --resume) RESUME=true; shift ;;
    --force) FORCE=true; shift ;;
    --dry-run) DRY_RUN=true; shift ;;
    -h|--help) usage; exit 0 ;;
    --workflow) die "--workflow is no longer supported. TITAN always runs the full pipeline." ;;
    --) shift; EXTRA_NF_ARGS+=("$@"); break ;;
    *) die "Unknown option: $1. Use --help." ;;
  esac
done

[[ -d "$PROJECT_DIR" ]] || die "Project directory does not exist: ${PROJECT_DIR}"
[[ -f "$PROJECT_DIR/main.nf" ]] || die "main.nf not found in project directory: ${PROJECT_DIR}"

require_file "previous assembly" "$PREVIOUS_ASSEMBLY"
require_file "new assembly" "$NEW_ASSEMBLY"
require_file "previous annotations" "$PREVIOUS_ANNOTATIONS"
require_file "RNA-seq samplesheet" "$RNASEQ_SAMPLESHEET"
require_dir "RNA-seq data directory" "$RNASEQ_DATA_DIR"
require_file "protein samplesheet" "$PROTEIN_SAMPLESHEET"
require_file "EGAPx parameter file" "$EGAPX_PARAMFILE"
[[ -n "$OUTPUT_DIR" ]] || die "output directory is required"
[[ -n "$PROFILE" ]] || die "profile is required"
[[ "$RUN_NAME" =~ ^[A-Za-z0-9_.-]+$ ]] || die "run name may only contain letters, numbers, dot, underscore and dash"

mkdir -p "$OUTPUT_DIR"
OUTPUT_DIR="$(cd "$OUTPUT_DIR" && pwd -P)"
WORK_DIR="${WORK_DIR:-$OUTPUT_DIR/work}"
EGAPX_CACHE_DIR="${EGAPX_CACHE_DIR:-$PROJECT_DIR/.egapx_cache}"
EGAPX_RUNNER_DIR="${EGAPX_RUNNER_DIR:-$PROJECT_DIR/.egapx_runner}"
REPORTS_DIR="$OUTPUT_DIR/nextflow_reports"
PREVIOUS_ASSEMBLY="$(abs_path "$PREVIOUS_ASSEMBLY")"
NEW_ASSEMBLY="$(abs_path "$NEW_ASSEMBLY")"
PREVIOUS_ANNOTATIONS="$(abs_path "$PREVIOUS_ANNOTATIONS")"
RNASEQ_SAMPLESHEET="$(abs_path "$RNASEQ_SAMPLESHEET")"
RNASEQ_DATA_DIR="$(abs_path "$RNASEQ_DATA_DIR")"
PROTEIN_SAMPLESHEET="$(abs_path "$PROTEIN_SAMPLESHEET")"
EGAPX_PARAMFILE="$(abs_path "$EGAPX_PARAMFILE")"
EGAPX_CACHE_DIR="$(abs_path "$EGAPX_CACHE_DIR")"
EGAPX_RUNNER_DIR="$(abs_path "$EGAPX_RUNNER_DIR")"

if [[ "$RESUME" == false && "$FORCE" == false ]] && find "$OUTPUT_DIR" -mindepth 1 -maxdepth 1 | grep -q .; then
  die "Output directory is not empty: ${OUTPUT_DIR}. Use --resume to continue or --force to start a new run there."
fi

mkdir -p "$WORK_DIR" "$REPORTS_DIR"

if [[ ",$PROFILE," == *",apptainer,"* ]]; then
  export TITAN_APPTAINER_CACHEDIR="${TITAN_APPTAINER_CACHEDIR:-$OUTPUT_DIR/apptainer-cache}"
  mkdir -p "$TITAN_APPTAINER_CACHEDIR"
fi

if [[ -n "$NEXTFLOW_MODULE" ]] && command -v module >/dev/null 2>&1; then
  module load "$NEXTFLOW_MODULE"
fi

command -v nextflow >/dev/null 2>&1 || die "nextflow is not available in PATH"
command -v python3 >/dev/null 2>&1 || die "python3 is not available in PATH"
command -v curl >/dev/null 2>&1 || die "curl is not available in PATH"
command -v tar >/dev/null 2>&1 || die "tar is not available in PATH"

cd "$PROJECT_DIR"

CONFIG_SNAPSHOT="$(mktemp)"
trap 'rm -f "$CONFIG_SNAPSHOT"' EXIT
nextflow config -profile "$PROFILE" > "$CONFIG_SNAPSHOT"

config_value() {
  local key="$1"
  awk -v key="$key" '
    $1 == key && $2 == "=" {
      value = $0
      sub("^[[:space:]]*" key "[[:space:]]*=[[:space:]]*", "", value)
      gsub(/^'\''|'\''$/, "", value)
      print value
      exit
    }
  ' "$CONFIG_SNAPSHOT"
}

python3 scripts/validate_container_pins.py >/dev/null
python3 scripts/validate_profiles.py >/dev/null
python3 scripts/validate_inputs.py \
  --project-dir "$PROJECT_DIR" \
  --new-assembly "$NEW_ASSEMBLY" \
  --previous-assembly "$PREVIOUS_ASSEMBLY" \
  --previous-annotations "$PREVIOUS_ANNOTATIONS" \
  --rnaseq-samplesheet "$RNASEQ_SAMPLESHEET" \
  --rnaseq-data-dir "$RNASEQ_DATA_DIR" \
  --protein-samplesheet "$PROTEIN_SAMPLESHEET" \
  --egapx-paramfile "$EGAPX_PARAMFILE" \
  --egapx-executor "$(config_value egapx_executor)" \
  --psiclass-vd "$(config_value PSICLASS_vd_option)" \
  --psiclass-c "$(config_value PSICLASS_c_option)" >/dev/null

prepare_egapx_cache() {
  local data_version="$1"
  local revision="$2"
  local runner_script="$EGAPX_RUNNER_DIR/ui/egapx.py"

  mkdir -p "$EGAPX_CACHE_DIR" "$EGAPX_RUNNER_DIR"
  if [[ ! -f "$runner_script" ]]; then
    find "$EGAPX_RUNNER_DIR" -mindepth 1 -maxdepth 1 -exec rm -rf {} +
    curl -fsSL "https://github.com/ncbi/egapx/archive/refs/tags/${revision}.tar.gz" \
      | tar -xz --strip-components=1 -C "$EGAPX_RUNNER_DIR"
  fi

  python3 "$runner_script" -dl -lc "$EGAPX_CACHE_DIR" -dv "$data_version"
}

if [[ "$PREPARE_EGAPX_CACHE" == true ]]; then
  prepare_egapx_cache "$(config_value egapx_data_version)" "$(config_value egapx_revision)"
fi

cmd=(
  nextflow run main.nf
  -profile "$PROFILE"
  -name "$RUN_NAME"
  -work-dir "$WORK_DIR"
  -with-report "$REPORTS_DIR/${RUN_NAME}.report.html"
  -with-timeline "$REPORTS_DIR/${RUN_NAME}.timeline.html"
  -with-trace "$REPORTS_DIR/${RUN_NAME}.trace.txt"
  -with-dag "$REPORTS_DIR/${RUN_NAME}.dag.html"
  -ansi-log false
  --output_dir "$OUTPUT_DIR"
  --previous_assembly "$PREVIOUS_ASSEMBLY"
  --new_assembly "$NEW_ASSEMBLY"
  --previous_annotations "$PREVIOUS_ANNOTATIONS"
  --RNAseq_samplesheet "$RNASEQ_SAMPLESHEET"
  --RNAseq_data_dir "$RNASEQ_DATA_DIR"
  --protein_samplesheet "$PROTEIN_SAMPLESHEET"
  --egapx_paramfile "$EGAPX_PARAMFILE"
  --egapx_runner_dir "$EGAPX_RUNNER_DIR"
  --egapx_local_cache_dir "$EGAPX_CACHE_DIR"
)

if [[ "$RESUME" == true ]]; then
  cmd+=(-resume)
fi

cmd+=("${EXTRA_NF_ARGS[@]}")

echo "Project directory: $PROJECT_DIR"
echo "Output directory:  $OUTPUT_DIR"
echo "Work directory:    $WORK_DIR"
echo "Reports directory: $REPORTS_DIR"
echo "Profile:           $PROFILE"
echo "Run name:          $RUN_NAME"
echo "EGAPx runner:      $EGAPX_RUNNER_DIR"
echo "EGAPx cache:       $EGAPX_CACHE_DIR"
if [[ "${TITAN_APPTAINER_CACHEDIR:-}" ]]; then
  echo "Apptainer cache:   $TITAN_APPTAINER_CACHEDIR"
fi
echo
echo "Command:"
quote_cmd "${cmd[@]}"

if [[ "$DRY_RUN" == true ]]; then
  exit 0
fi

exec "${cmd[@]}"
