#!/usr/bin/env bash
#SBATCH --job-name=mikado_test
#SBATCH --nodelist=node005
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=%x-%j.log
# Orchestrator for one arm of the Mikado BRAKER3-input test (see mikado_input_test.nf).
# Like launch_TITAN_serveur_colmar.sh, this job only runs Nextflow; every task is its own sbatch job.
#
#   sbatch -x calcul scripts/mikado_input_test/launch_arm.sh raw       # control: production inputs
#   sbatch -x calcul scripts/mikado_input_test/launch_arm.sh tsebra    # test: braker.gff3 + genemark_supported.gtf
# Mikado tasks get 24 CPUs (scripts/mikado_input_test/resources.config).
#
# Add -resume automatically when the arm was already started.  Outputs: data/mikado_input_test/<arm>/
set -Eeuo pipefail
ARM=${1:?usage: launch_arm.sh raw|tsebra}
[[ $ARM == raw || $ARM == tsebra ]] || { echo "arm must be raw or tsebra" >&2; exit 2; }

PROJECT_DIR="${TITAN_PROJECT_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd -P)}}"
CONFIG_FILE="${TITAN_CONFIG_FILE:-$PROJECT_DIR/data/slurm_apptainer.config}"
OUT="$PROJECT_DIR/data/mikado_input_test/$ARM"
WORK="$PROJECT_DIR/data/mikado_input_test/work_$ARM"
cd "$PROJECT_DIR"
[[ -f mikado_input_test.nf && -f $CONFIG_FILE ]] || { echo "run from the TITAN project root" >&2; exit 2; }

# `module load` is unreliable in a non-login batch shell (on node005 it leaves Java 1.8 active and
# the python module calls a missing /usr/bin/scl), so call Lmod directly; python is not needed here.
if [[ -n ${LMOD_CMD:-} ]]; then
  eval "$("$LMOD_CMD" bash load Java/17.0.13 nextflow/24.04.3 apptainer/1.4.0-rc.2)"
fi
command -v nextflow >/dev/null && command -v apptainer >/dev/null || { echo "nextflow/apptainer not available" >&2; exit 2; }
java -version 2>&1 | grep -q '"1[1-9]\|"2[0-9]' || { echo "Java >= 11 required" >&2; exit 2; }
mkdir -p "$OUT/05_run_info/nextflow_reports" "$WORK" "$PROJECT_DIR/.tmp" "$PROJECT_DIR/.apptainer-tmp"
export NXF_HOME="${NXF_HOME:-$PROJECT_DIR/.nextflow-home}"
export TMPDIR="$PROJECT_DIR/.tmp"
export APPTAINER_CACHEDIR="$PROJECT_DIR/.apptainer-cache" SINGULARITY_CACHEDIR="$PROJECT_DIR/.apptainer-cache"
export APPTAINER_TMPDIR="$PROJECT_DIR/.apptainer-tmp" SINGULARITY_TMPDIR="$PROJECT_DIR/.apptainer-tmp"
export PYTHONNOUSERSITE=1 DEBUGINFOD_URLS=/dev/null

resume=()
[[ -d $WORK/.. && -n $(ls -A "$WORK" 2>/dev/null) ]] && resume=(-resume)

nextflow -c "$CONFIG_FILE" -c "$PROJECT_DIR/scripts/mikado_input_test/resources.config" run mikado_input_test.nf \
  -profile slurm,apptainer \
  -name "mikado_test_${ARM}_$(date +%Y%m%d_%H%M%S)" \
  -work-dir "$WORK" \
  -ansi-log false \
  "${resume[@]}" \
  --braker_mode "$ARM" \
  --output_dir "$OUT"
