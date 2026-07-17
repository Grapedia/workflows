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
  --prepare-eggnog-data          Download the eggNOG-mapper database before launching TITAN
                                 Default: TITAN_PREPARE_EGGNOG_DATA or false
  --eggnog-data-dir PATH         Persistent eggNOG-mapper database directory
                                 Default: TITAN_EGGNOG_DATA_DIR or <project-dir>/.eggnog_data
  --enable-eggnog-mapper         Pass --run_eggnog_mapper true and --eggnog_data_dir to Nextflow
                                 Default: TITAN_RUN_EGGNOG_MAPPER or false
  --prepare-helixer-model        Download the Helixer lineage model before launching TITAN
                                 Default: TITAN_PREPARE_HELIXER_MODEL or false
  --helixer-model-dir PATH       Persistent Helixer model directory
                                 Default: TITAN_HELIXER_MODEL_DIR or <project-dir>/.helixer_models
  --helixer-lineage NAME         Helixer lineage to fetch/use: vertebrate, land_plant, fungi or invertebrate
                                 Default: TITAN_HELIXER_LINEAGE or land_plant
  --enable-helixer               Pass --run_helixer true, --helixer_model_dir and --helixer_model to Nextflow
                                 Default: TITAN_RUN_HELIXER or false
  --enable-helixer-gpu           Pass --helixer_use_gpu true to Nextflow (requires a GPU visible on the node)
                                 Default: TITAN_HELIXER_USE_GPU or false
  --prepare-interproscan-data    Download the InterProScan member database data before launching TITAN
                                 Default: TITAN_PREPARE_INTERPROSCAN_DATA or false
  --interproscan-data-dir PATH   Persistent InterProScan member database directory
                                 Default: TITAN_INTERPROSCAN_DATA_DIR or <project-dir>/.interproscan_data
  --enable-interproscan          Pass --run_interproscan true and --interproscan_data_dir to Nextflow
                                 Default: TITAN_RUN_INTERPROSCAN or false
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
PREPARE_EGGNOG_DATA="${TITAN_PREPARE_EGGNOG_DATA:-false}"
EGGNOG_DATA_DIR="${TITAN_EGGNOG_DATA_DIR:-}"
RUN_EGGNOG_MAPPER="${TITAN_RUN_EGGNOG_MAPPER:-false}"
PREPARE_HELIXER_MODEL="${TITAN_PREPARE_HELIXER_MODEL:-false}"
HELIXER_MODEL_DIR="${TITAN_HELIXER_MODEL_DIR:-}"
HELIXER_LINEAGE="${TITAN_HELIXER_LINEAGE:-land_plant}"
RUN_HELIXER="${TITAN_RUN_HELIXER:-false}"
HELIXER_USE_GPU="${TITAN_HELIXER_USE_GPU:-false}"
PREPARE_INTERPROSCAN_DATA="${TITAN_PREPARE_INTERPROSCAN_DATA:-false}"
INTERPROSCAN_DATA_DIR="${TITAN_INTERPROSCAN_DATA_DIR:-}"
RUN_INTERPROSCAN="${TITAN_RUN_INTERPROSCAN:-false}"
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
    --prepare-eggnog-data) PREPARE_EGGNOG_DATA=true; shift ;;
    --eggnog-data-dir) EGGNOG_DATA_DIR="${2:-}"; shift 2 ;;
    --enable-eggnog-mapper) RUN_EGGNOG_MAPPER=true; shift ;;
    --prepare-helixer-model) PREPARE_HELIXER_MODEL=true; shift ;;
    --helixer-model-dir) HELIXER_MODEL_DIR="${2:-}"; shift 2 ;;
    --helixer-lineage) HELIXER_LINEAGE="${2:-}"; shift 2 ;;
    --enable-helixer) RUN_HELIXER=true; shift ;;
    --enable-helixer-gpu) HELIXER_USE_GPU=true; shift ;;
    --prepare-interproscan-data) PREPARE_INTERPROSCAN_DATA=true; shift ;;
    --interproscan-data-dir) INTERPROSCAN_DATA_DIR="${2:-}"; shift 2 ;;
    --enable-interproscan) RUN_INTERPROSCAN=true; shift ;;
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
EGGNOG_DATA_DIR="${EGGNOG_DATA_DIR:-$PROJECT_DIR/.eggnog_data}"
HELIXER_MODEL_DIR="${HELIXER_MODEL_DIR:-$PROJECT_DIR/.helixer_models}"
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
mkdir -p "$EGGNOG_DATA_DIR"
EGGNOG_DATA_DIR="$(abs_path "$EGGNOG_DATA_DIR")"
mkdir -p "$HELIXER_MODEL_DIR"
HELIXER_MODEL_DIR="$(abs_path "$HELIXER_MODEL_DIR")"
INTERPROSCAN_DATA_DIR="${INTERPROSCAN_DATA_DIR:-$PROJECT_DIR/.interproscan_data}"
mkdir -p "$INTERPROSCAN_DATA_DIR"
INTERPROSCAN_DATA_DIR="$(abs_path "$INTERPROSCAN_DATA_DIR")"

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

prepare_eggnog_data() {
  "$PROJECT_DIR/scripts/download_eggnog_data.sh" --data-dir "$EGGNOG_DATA_DIR"
}

if [[ "$PREPARE_EGGNOG_DATA" == true ]]; then
  prepare_eggnog_data
fi

prepare_helixer_model() {
  local container="$1"
  "$PROJECT_DIR/scripts/download_helixer_model.sh" \
    --model-dir "$HELIXER_MODEL_DIR" \
    --container "$container" \
    --lineage "$HELIXER_LINEAGE"
}

if [[ "$PREPARE_HELIXER_MODEL" == true ]]; then
  prepare_helixer_model "$(config_value container_helixer)"
fi

prepare_interproscan_data() {
  "$PROJECT_DIR/scripts/download_interproscan_data.sh" --data-dir "$INTERPROSCAN_DATA_DIR"
}

if [[ "$PREPARE_INTERPROSCAN_DATA" == true ]]; then
  prepare_interproscan_data
fi

cmd=(
  nextflow run main.nf
  -profile "$PROFILE"
  -name "$RUN_NAME"
  -work-dir "$WORK_DIR"
  # -with-report/-with-timeline/-with-trace are intentionally omitted: they
  # require `ps` inside every task's execution context to sample runtime
  # metrics, and none of TITAN's pinned container images ship procps (and
  # Apptainer does not leak host binaries into the container filesystem), so
  # every single task fails immediately with "Command 'ps' required by
  # nextflow to collect task metrics cannot be found" when any of them are
  # enabled. -with-dag does not need runtime metrics and is safe to keep.
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

if [[ "$RUN_EGGNOG_MAPPER" == true ]]; then
  cmd+=(--run_eggnog_mapper true --eggnog_data_dir "$EGGNOG_DATA_DIR")
fi

if [[ "$RUN_HELIXER" == true ]]; then
  cmd+=(--run_helixer true --helixer_model_dir "$HELIXER_MODEL_DIR" --helixer_model "$HELIXER_LINEAGE")
  if [[ "$HELIXER_USE_GPU" == true ]]; then
    cmd+=(--helixer_use_gpu true)
  fi
fi

if [[ "$RUN_INTERPROSCAN" == true ]]; then
  cmd+=(--run_interproscan true --interproscan_data_dir "$INTERPROSCAN_DATA_DIR")
fi

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
if [[ "$RUN_EGGNOG_MAPPER" == true ]]; then
  echo "eggNOG-mapper:     enabled, data dir $EGGNOG_DATA_DIR"
fi
if [[ "$RUN_HELIXER" == true ]]; then
  echo "Helixer:           enabled, lineage $HELIXER_LINEAGE, model dir $HELIXER_MODEL_DIR, gpu $HELIXER_USE_GPU"
fi
if [[ "$RUN_INTERPROSCAN" == true ]]; then
  echo "InterProScan:      enabled, data dir $INTERPROSCAN_DATA_DIR"
fi
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
