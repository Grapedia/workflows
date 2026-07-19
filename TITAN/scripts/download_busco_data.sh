#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
  cat <<'EOF'
Usage:
  scripts/download_busco_data.sh --data-dir PATH --container IMAGE --lineage LINEAGE [--runtime apptainer|singularity|docker]

Downloads a BUSCO lineage dataset into --data-dir by running `busco --download`
inside the same pinned BUSCO container image used by TITAN.

Pass the resulting directory as --busco_data_dir to TITAN, together with
--run_busco true and --busco_lineage <LINEAGE>. TITAN later runs BUSCO with
--offline --download_path <data-dir>, so the production run does not need to
download reference data.

Options:
  --data-dir PATH  Target directory for BUSCO downloaded lineages. Created if missing.
  --container IMAGE Pinned BUSCO image reference, e.g. the value of
                    params.container_busco.
  --lineage NAME   BUSCO lineage dataset to fetch, e.g. eudicotyledons_odb12.2.
  --runtime NAME   Container runtime to use: apptainer, singularity or docker.
                    Default: first of apptainer, singularity, docker found in PATH.
  -h, --help       Show this help
EOF
}

die() {
  echo "ERROR: $*" >&2
  exit 1
}

DATA_DIR=""
CONTAINER=""
LINEAGE=""
RUNTIME=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --data-dir) DATA_DIR="${2:-}"; shift 2 ;;
    --container) CONTAINER="${2:-}"; shift 2 ;;
    --lineage) LINEAGE="${2:-}"; shift 2 ;;
    --runtime) RUNTIME="${2:-}"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown option: $1. Use --help." ;;
  esac
done

[[ -n "$DATA_DIR" ]] || die "--data-dir is required"
[[ -n "$CONTAINER" ]] || die "--container is required"
[[ -n "$LINEAGE" ]] || die "--lineage is required"

mkdir -p "$DATA_DIR"
DATA_DIR="$(cd "$DATA_DIR" && pwd -P)"

lineage_present() {
  [[ -d "${DATA_DIR}/lineages/${LINEAGE}" ]] || [[ -d "${DATA_DIR}/${LINEAGE}" ]]
}

if lineage_present; then
  echo "BUSCO lineage already present, skipping: ${LINEAGE} in ${DATA_DIR}"
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

echo "Downloading BUSCO lineage '${LINEAGE}' into ${DATA_DIR} using ${RUNTIME} (${CONTAINER})"

case "$RUNTIME" in
  apptainer|singularity)
    "$RUNTIME" exec --bind "${DATA_DIR}:${DATA_DIR}" "docker://${CONTAINER}" \
      busco --download "${LINEAGE}" --download_path "${DATA_DIR}"
    ;;
  docker)
    docker run --rm -v "${DATA_DIR}:${DATA_DIR}" "${CONTAINER}" \
      busco --download "${LINEAGE}" --download_path "${DATA_DIR}"
    ;;
  *)
    die "unsupported runtime: ${RUNTIME}"
    ;;
esac

lineage_present || die "download appears incomplete: ${LINEAGE} not found under ${DATA_DIR} or ${DATA_DIR}/lineages"

echo "BUSCO lineage '${LINEAGE}' ready under ${DATA_DIR}"
