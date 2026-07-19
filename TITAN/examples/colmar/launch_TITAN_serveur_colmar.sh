#!/usr/bin/env bash
#SBATCH --job-name=TITAN_colmar
#SBATCH --nodelist=calcul,node001,node003,node004,node005
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=%x-%j.log
# This script is itself only the Nextflow orchestrator: it stays alive for the
# whole pipeline (potentially several days across EGAPx/BRAKER3/InterProScan)
# and submits every actual TITAN task as its own separate sbatch job via
# Nextflow's slurm executor (see clusterOptions/--nodelist in
# data/slurm_apptainer.config, which restricts those task jobs the same way).
# No --time is set here on purpose: an incorrect guess could make sbatch
# reject the submission outright if it exceeds your partition's default. If
# your site's default walltime is too short for the full run, submit with
# e.g. `sbatch --time=14-00:00:00 launch_TITAN_serveur_colmar.sh` instead.
set -Eeuo pipefail


# Under sbatch, Slurm copies this script into a per-job spool directory before
# running it, so ${BASH_SOURCE[0]} no longer points at the real script
# location; SLURM_SUBMIT_DIR (the cwd at `sbatch` time) is the reliable path
# in that case. Fall back to BASH_SOURCE-based detection for direct/local runs.
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
  SCRIPT_DIR="$SLURM_SUBMIT_DIR"
else
  SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
fi
PROJECT_DIR="${TITAN_PROJECT_DIR:-$SCRIPT_DIR}"
RUN_NAME="${TITAN_RUN_NAME:-PN40024_T2T_prod_$(date +%Y%m%d_%H%M%S)}"
OUTPUT_DIR="${TITAN_OUTPUT_DIR:-$PROJECT_DIR/data/titan_prod_out}"
WORK_DIR="${TITAN_WORK_DIR:-$PROJECT_DIR/data/work}"
CONFIG_FILE="${TITAN_CONFIG_FILE:-$PROJECT_DIR/data/slurm_apptainer.config}"
PROFILE="${TITAN_PROFILE:-slurm,apptainer}"
# Plain `-resume` (no argument) resumes from whatever session is most recent
# in .nextflow/history - including unrelated one-off `nextflow run` commands
# (e.g. a `-stub-run`/`-preview` sanity check) run from this same project
# directory after the last real production attempt. That silently detaches
# the resume chain from the actual production cache and reruns everything
# from scratch. Set TITAN_RESUME_ID to the session UUID of the last real
# production run (see the 6th column of `.nextflow/history`, or `nextflow
# log`) to pin resume to it explicitly instead of trusting "most recent".
RESUME_ID="${TITAN_RESUME_ID:-}"
PREPARE_EGAPX_CACHE="${TITAN_PREPARE_EGAPX_CACHE:-false}"
EGAPX_CACHE_DIR="${TITAN_EGAPX_CACHE_DIR:-$PROJECT_DIR/.egapx_cache}"
EGAPX_RUNNER_DIR="${TITAN_EGAPX_RUNNER_DIR:-$PROJECT_DIR/.egapx_runner}"
EGAPX_REVISION="${TITAN_EGAPX_REVISION:-v0.5.2}"
EGAPX_DATA_VERSION="${TITAN_EGAPX_DATA_VERSION:-egapxsupportdata_20251017}"
PREPARE_EGGNOG_DATA="${TITAN_PREPARE_EGGNOG_DATA:-false}"
EGGNOG_DATA_DIR="${TITAN_EGGNOG_DATA_DIR:-$PROJECT_DIR/.eggnog_data}"
PREPARE_HELIXER_MODEL="${TITAN_PREPARE_HELIXER_MODEL:-false}"
HELIXER_MODEL_DIR="${TITAN_HELIXER_MODEL_DIR:-$PROJECT_DIR/.helixer_models}"
HELIXER_LINEAGE="${TITAN_HELIXER_LINEAGE:-land_plant}"
PREPARE_INTERPROSCAN_DATA="${TITAN_PREPARE_INTERPROSCAN_DATA:-false}"
INTERPROSCAN_DATA_DIR="${TITAN_INTERPROSCAN_DATA_DIR:-$PROJECT_DIR/.interproscan_data}"
PREPARE_BUSCO_DATA="${TITAN_PREPARE_BUSCO_DATA:-false}"
BUSCO_DATA_DIR="${TITAN_BUSCO_DATA_DIR:-$PROJECT_DIR/.busco_data}"
BUSCO_LINEAGE="${TITAN_BUSCO_LINEAGE:-eudicotyledons_odb12.2}"
PREPARE_RFAM_DATA="${TITAN_PREPARE_RFAM_DATA:-false}"
RFAM_DATA_DIR="${TITAN_RFAM_DATA_DIR:-$PROJECT_DIR/.rfam_data}"
PREPARE_OMARK_DATA="${TITAN_PREPARE_OMARK_DATA:-false}"
OMARK_DATA_DIR="${TITAN_OMARK_DATA_DIR:-$PROJECT_DIR/.omark_data}"

usage() {
  cat <<'EOF'
Usage:
  ./launch_TITAN_serveur_colmar.sh [--prepare-egapx-cache] [--prepare-eggnog-data] [--prepare-helixer-model] [--prepare-interproscan-data] [--prepare-busco-data] [--prepare-rfam-data] [--prepare-omark-data] [--dry-run] [-- no_more_nextflow_args]

Environment overrides:
  TITAN_OUTPUT_DIR, TITAN_WORK_DIR, TITAN_RUN_NAME, TITAN_CONFIG_FILE,
  TITAN_EGAPX_CACHE_DIR, TITAN_EGAPX_RUNNER_DIR, TITAN_PREPARE_EGAPX_CACHE,
  TITAN_EGGNOG_DATA_DIR, TITAN_PREPARE_EGGNOG_DATA,
  TITAN_HELIXER_MODEL_DIR, TITAN_HELIXER_LINEAGE, TITAN_PREPARE_HELIXER_MODEL,
  TITAN_INTERPROSCAN_DATA_DIR, TITAN_PREPARE_INTERPROSCAN_DATA,
  TITAN_BUSCO_DATA_DIR, TITAN_BUSCO_LINEAGE, TITAN_PREPARE_BUSCO_DATA,
  TITAN_RFAM_DATA_DIR, TITAN_PREPARE_RFAM_DATA,
  TITAN_OMARK_DATA_DIR, TITAN_PREPARE_OMARK_DATA.

Default production inputs are defined in data/slurm_apptainer.config. To run
eggNOG-mapper in production, set run_eggnog_mapper = true in that config file
and pass --prepare-eggnog-data at least once to populate TITAN_EGGNOG_DATA_DIR
(default <project-dir>/.eggnog_data).

To run Helixer in production, set run_helixer = true (and optionally
helixer_use_gpu = true, only if a GPU is available on the node) in that
config file, and pass --prepare-helixer-model at least once to populate
TITAN_HELIXER_MODEL_DIR (default <project-dir>/.helixer_models).

To run InterProScan in production, set run_interproscan = true in that config
file and pass --prepare-interproscan-data at least once to populate
TITAN_INTERPROSCAN_DATA_DIR (default <project-dir>/.interproscan_data). This
is a large download (several GB compressed); it only needs to run once.

To run Infernal/Rfam in production, set run_rfam = true in that config file
and pass --prepare-rfam-data at least once to populate TITAN_RFAM_DATA_DIR
(default <project-dir>/.rfam_data) with the indexed Rfam.cm/.clanin library.

To run BUSCO in production, set run_busco = true in that config file and pass
--prepare-busco-data at least once to populate TITAN_BUSCO_DATA_DIR (default
<project-dir>/.busco_data) with the selected TITAN_BUSCO_LINEAGE dataset.

To run OMArk in production, set run_omark = true in that config file and
pass --prepare-omark-data at least once to populate TITAN_OMARK_DATA_DIR
(default <project-dir>/.omark_data) with the OMAmer LUCA.h5 database. This is
a large download (~9.5 GB); it only needs to run once.
EOF
}

die() {
  echo "ERROR: $*" >&2
  exit 1
}

quote_cmd() {
  printf '%q ' "$@"
  printf '\n'
}

DRY_RUN=false
EXTRA_NF_ARGS=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --prepare-egapx-cache) PREPARE_EGAPX_CACHE=true; shift ;;
    --prepare-eggnog-data) PREPARE_EGGNOG_DATA=true; shift ;;
    --prepare-helixer-model) PREPARE_HELIXER_MODEL=true; shift ;;
    --prepare-interproscan-data) PREPARE_INTERPROSCAN_DATA=true; shift ;;
    --prepare-busco-data) PREPARE_BUSCO_DATA=true; shift ;;
    --prepare-rfam-data) PREPARE_RFAM_DATA=true; shift ;;
    --prepare-omark-data) PREPARE_OMARK_DATA=true; shift ;;
    --dry-run) DRY_RUN=true; shift ;;
    -h|--help) usage; exit 0 ;;
    --) shift; EXTRA_NF_ARGS+=("$@"); break ;;
    *) die "Unknown option: $1. Use --help." ;;
  esac
done

cd "$PROJECT_DIR"

[[ -f "$CONFIG_FILE" ]] || die "config file not found: $CONFIG_FILE"
[[ -f main.nf ]] || die "main.nf not found in $PROJECT_DIR"

if command -v module >/dev/null 2>&1; then
  module load Java/17.0.13
  module load nextflow/24.04.3
  module load apptainer/1.4.0-rc.2
  module load python/3.12
fi

command -v nextflow >/dev/null 2>&1 || die "nextflow is not available in PATH"
command -v apptainer >/dev/null 2>&1 || die "apptainer is not available in PATH"
command -v python3 >/dev/null 2>&1 || die "python3 is not available in PATH"
command -v curl >/dev/null 2>&1 || die "curl is not available in PATH"
command -v tar >/dev/null 2>&1 || die "tar is not available in PATH"

mkdir -p \
  "$OUTPUT_DIR/nextflow_reports" \
  "$WORK_DIR" \
  "$EGAPX_CACHE_DIR" \
  "$EGAPX_RUNNER_DIR" \
  "$EGGNOG_DATA_DIR" \
  "$HELIXER_MODEL_DIR" \
  "$INTERPROSCAN_DATA_DIR" \
  "$BUSCO_DATA_DIR" \
  "$RFAM_DATA_DIR" \
  "$OMARK_DATA_DIR" \
  "$PROJECT_DIR/.apptainer-cache" \
  "$PROJECT_DIR/.apptainer-tmp" \
  "$PROJECT_DIR/.nextflow-home" \
  "$PROJECT_DIR/.tmp"

export NXF_HOME="${NXF_HOME:-$PROJECT_DIR/.nextflow-home}"
export TMPDIR="${TMPDIR:-$PROJECT_DIR/.tmp}"
export TITAN_APPTAINER_CACHEDIR="${TITAN_APPTAINER_CACHEDIR:-$PROJECT_DIR/.apptainer-cache}"
export APPTAINER_CACHEDIR="${APPTAINER_CACHEDIR:-$TITAN_APPTAINER_CACHEDIR}"
export SINGULARITY_CACHEDIR="${SINGULARITY_CACHEDIR:-$TITAN_APPTAINER_CACHEDIR}"
export APPTAINER_TMPDIR="${APPTAINER_TMPDIR:-$PROJECT_DIR/.apptainer-tmp}"
export SINGULARITY_TMPDIR="${SINGULARITY_TMPDIR:-$PROJECT_DIR/.apptainer-tmp}"
export PYTHONNOUSERSITE=1
export DEBUGINFOD_URLS=/dev/null

CONFIG_SNAPSHOT="$(mktemp)"
trap 'rm -f "$CONFIG_SNAPSHOT"' EXIT
nextflow -c "$CONFIG_FILE" config -profile "$PROFILE" > "$CONFIG_SNAPSHOT"

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

prepare_egapx_cache() {
  local runner_script="$EGAPX_RUNNER_DIR/ui/egapx.py"

  if [[ ! -f "$runner_script" ]]; then
    find "$EGAPX_RUNNER_DIR" -mindepth 1 -maxdepth 1 -exec rm -rf {} +
    curl -fsSL "https://github.com/ncbi/egapx/archive/refs/tags/${EGAPX_REVISION}.tar.gz" \
      | tar -xz --strip-components=1 -C "$EGAPX_RUNNER_DIR"
  fi

  python3 "$runner_script" -dl -lc "$EGAPX_CACHE_DIR" -dv "$EGAPX_DATA_VERSION" -dn "29760"
}

python3 scripts/validate_container_pins.py >/dev/null
python3 scripts/validate_profiles.py >/dev/null
python3 scripts/validate_inputs.py \
  --project-dir "$PROJECT_DIR" \
  --new-assembly "$(config_value new_assembly)" \
  --previous-assembly "$(config_value previous_assembly)" \
  --previous-annotations "$(config_value previous_annotations)" \
  --rnaseq-samplesheet "$(config_value RNAseq_samplesheet)" \
  --rnaseq-data-dir "$(config_value RNAseq_data_dir)" \
  --protein-samplesheet "$(config_value protein_samplesheet)" \
  --egapx-paramfile "$(config_value egapx_paramfile)" \
  --egapx-executor "$(config_value egapx_executor)" \
  --psiclass-vd "$(config_value PSICLASS_vd_option)" \
  --psiclass-c "$(config_value PSICLASS_c_option)" >/dev/null

if [[ "$PREPARE_EGAPX_CACHE" == true ]]; then
  prepare_egapx_cache
fi

prepare_eggnog_data() {
  "$PROJECT_DIR/scripts/download_eggnog_data.sh" --data-dir "$EGGNOG_DATA_DIR"
}

if [[ "$PREPARE_EGGNOG_DATA" == true ]]; then
  prepare_eggnog_data
fi

prepare_helixer_model() {
  "$PROJECT_DIR/scripts/download_helixer_model.sh" \
    --model-dir "$HELIXER_MODEL_DIR" \
    --container "$(config_value container_helixer)" \
    --lineage "$HELIXER_LINEAGE"
}

if [[ "$PREPARE_HELIXER_MODEL" == true ]]; then
  prepare_helixer_model
fi

prepare_interproscan_data() {
  "$PROJECT_DIR/scripts/download_interproscan_data.sh" --data-dir "$INTERPROSCAN_DATA_DIR"
}

if [[ "$PREPARE_INTERPROSCAN_DATA" == true ]]; then
  prepare_interproscan_data
fi

prepare_busco_data() {
  "$PROJECT_DIR/scripts/download_busco_data.sh" \
    --data-dir "$BUSCO_DATA_DIR" \
    --container "$(config_value container_busco)" \
    --lineage "$BUSCO_LINEAGE"
}

if [[ "$PREPARE_BUSCO_DATA" == true ]]; then
  prepare_busco_data
fi

prepare_rfam_data() {
  "$PROJECT_DIR/scripts/download_rfam_data.sh" \
    --data-dir "$RFAM_DATA_DIR" \
    --container "$(config_value container_infernal)"
}

if [[ "$PREPARE_RFAM_DATA" == true ]]; then
  prepare_rfam_data
fi

prepare_omark_data() {
  "$PROJECT_DIR/scripts/download_omark_data.sh" --data-dir "$OMARK_DATA_DIR"
}

if [[ "$PREPARE_OMARK_DATA" == true ]]; then
  prepare_omark_data
fi

resume_args=(-resume)
if [[ -n "$RESUME_ID" ]]; then
  resume_args=(-resume "$RESUME_ID")
fi

cmd=(
  nextflow -c "$CONFIG_FILE" run main.nf
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
  -with-dag "$OUTPUT_DIR/nextflow_reports/${RUN_NAME}.dag.html"
  -ansi-log false
  "${resume_args[@]}"
  --egapx_runner_dir "$EGAPX_RUNNER_DIR"
  --egapx_local_cache_dir "$EGAPX_CACHE_DIR"
  --eggnog_data_dir "$EGGNOG_DATA_DIR"
  --helixer_model_dir "$HELIXER_MODEL_DIR"
  --helixer_model "$HELIXER_LINEAGE"
  --interproscan_data_dir "$INTERPROSCAN_DATA_DIR"
  --busco_data_dir "$BUSCO_DATA_DIR"
  --busco_lineage "$BUSCO_LINEAGE"
  --rfam_data_dir "$RFAM_DATA_DIR"
  --omark_data_dir "$OMARK_DATA_DIR"
)

cmd+=("${EXTRA_NF_ARGS[@]}")

echo "Project directory: $PROJECT_DIR"
echo "Config file:       $CONFIG_FILE"
echo "Output directory:  $OUTPUT_DIR"
echo "Work directory:    $WORK_DIR"
echo "Profile:           $PROFILE"
echo "Run name:          $RUN_NAME"
echo "Apptainer cache:   $APPTAINER_CACHEDIR"
echo "EGAPx cache:       $EGAPX_CACHE_DIR"
echo "eggNOG data dir:   $EGGNOG_DATA_DIR"
echo "Helixer model dir: $HELIXER_MODEL_DIR (lineage: $HELIXER_LINEAGE)"
echo "InterProScan data dir: $INTERPROSCAN_DATA_DIR"
echo "BUSCO data dir:    $BUSCO_DATA_DIR (lineage: $BUSCO_LINEAGE)"
echo "Rfam data dir:     $RFAM_DATA_DIR"
echo "OMArk data dir:    $OMARK_DATA_DIR"
echo "Resume ID:         ${RESUME_ID:-<latest Nextflow session>}"
echo
echo "Command:"
quote_cmd "${cmd[@]}"

if [[ "$DRY_RUN" == true ]]; then
  exit 0
fi

exec "${cmd[@]}"
