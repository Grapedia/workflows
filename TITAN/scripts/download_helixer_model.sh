#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
  cat <<'EOF'
Usage:
  scripts/download_helixer_model.sh --model-dir PATH --container IMAGE [--lineage LINEAGE] [--runtime apptainer|singularity|docker]

Downloads a pre-trained Helixer lineage model into --model-dir using the
pinned Helixer container's own fetch_helixer_models.py script. Helixer.py
itself does NOT auto-download models and will fail fast if the model is
missing; this script fetches it ahead of time so production runs stay
offline-reproducible.

Pass the resulting directory as --helixer_model_dir to TITAN, together with
--run_helixer true and --helixer_model <lineage>.

Options:
  --model-dir PATH  Target directory for Helixer models. Created if missing.
  --container IMAGE Pinned Helixer image reference, e.g. the value of
                    params.container_helixer.
  --lineage NAME    Helixer lineage to fetch: vertebrate, land_plant, fungi or
                    invertebrate. Default: land_plant.
  --runtime NAME    Container runtime to use: apptainer, singularity or docker.
                    Default: first of apptainer, singularity, docker found in PATH.
  -h, --help        Show this help
EOF
}

die() {
  echo "ERROR: $*" >&2
  exit 1
}

MODEL_DIR=""
CONTAINER=""
LINEAGE="land_plant"
RUNTIME=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --model-dir) MODEL_DIR="${2:-}"; shift 2 ;;
    --container) CONTAINER="${2:-}"; shift 2 ;;
    --lineage) LINEAGE="${2:-}"; shift 2 ;;
    --runtime) RUNTIME="${2:-}"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown option: $1. Use --help." ;;
  esac
done

[[ -n "$MODEL_DIR" ]] || die "--model-dir is required"
[[ -n "$CONTAINER" ]] || die "--container is required"

model_present() {
  find "${MODEL_DIR}/${LINEAGE}" -maxdepth 1 -name '*.h5' -size +0c 2>/dev/null | grep -q .
}

mkdir -p "$MODEL_DIR"
MODEL_DIR="$(cd "$MODEL_DIR" && pwd -P)"

if model_present; then
  echo "Already present, skipping: ${LINEAGE} model in ${MODEL_DIR}"
  exit 0
fi

if [[ -z "$RUNTIME" ]]; then
  if command -v apptainer >/dev/null 2>&1; then
    RUNTIME=apptainer
  elif command -v singularity >/dev/null 2>&1; then
    RUNTIME=singularity
  elif command -v docker >/dev/null 2>&1; then
    RUNTIME=docker
  else
    die "no container runtime found in PATH (apptainer, singularity or docker); load one or pass --runtime"
  fi
fi

echo "Fetching Helixer '${LINEAGE}' model into ${MODEL_DIR} using ${RUNTIME} (${CONTAINER})"

case "$RUNTIME" in
  apptainer|singularity)
    "$RUNTIME" exec --bind "${MODEL_DIR}:${MODEL_DIR}" "docker://${CONTAINER}" \
      fetch_helixer_models.py --lineage "${LINEAGE}" --custom-path "${MODEL_DIR}"
    ;;
  docker)
    docker run --rm -v "${MODEL_DIR}:${MODEL_DIR}" "${CONTAINER}" \
      fetch_helixer_models.py --lineage "${LINEAGE}" --custom-path "${MODEL_DIR}"
    ;;
  *)
    die "unsupported runtime: ${RUNTIME}"
    ;;
esac

model_present || die "download appears incomplete: no .h5 model found in ${MODEL_DIR}/${LINEAGE}"

echo "Helixer '${LINEAGE}' model ready at ${MODEL_DIR}/${LINEAGE}"
