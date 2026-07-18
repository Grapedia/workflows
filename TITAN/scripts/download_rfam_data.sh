#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
  cat <<'EOF'
Usage:
  scripts/download_rfam_data.sh --data-dir PATH --container IMAGE [--base-url URL] [--runtime apptainer|singularity|docker]

Downloads the Rfam covariance model library (Rfam.cm + Rfam.clanin) into
--data-dir with plain curl/gunzip, then indexes it with `cmpress` (run inside
the pinned Infernal container, since cmpress itself is an Infernal binary).

Pass the resulting directory as --rfam_data_dir to TITAN, together with
--run_rfam true.

Options:
  --data-dir PATH  Target directory for the Rfam.cm/.clanin/.i1* files. Created if missing.
  --container IMAGE Pinned Infernal image reference, e.g. the value of
                    params.container_infernal.
  --base-url URL   Base URL serving Rfam.cm.gz and Rfam.clanin.
                   Default: https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT
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
BASE_URL="https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT"
RUNTIME=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --data-dir) DATA_DIR="${2:-}"; shift 2 ;;
    --container) CONTAINER="${2:-}"; shift 2 ;;
    --base-url) BASE_URL="${2:-}"; shift 2 ;;
    --runtime) RUNTIME="${2:-}"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown option: $1. Use --help." ;;
  esac
done

[[ -n "$DATA_DIR" ]] || die "--data-dir is required"
[[ -n "$CONTAINER" ]] || die "--container is required"

command -v curl >/dev/null 2>&1 || die "curl is not available in PATH"
command -v gunzip >/dev/null 2>&1 || die "gunzip is not available in PATH"

mkdir -p "$DATA_DIR"
DATA_DIR="$(cd "$DATA_DIR" && pwd -P)"

index_present() {
  [[ -s "${DATA_DIR}/Rfam.cm.i1m" ]]
}

if index_present; then
  echo "Rfam data already present and indexed, skipping: ${DATA_DIR}"
  exit 0
fi

echo "Downloading Rfam library into ${DATA_DIR} from ${BASE_URL}"

if [[ ! -s "${DATA_DIR}/Rfam.cm" ]]; then
  curl -fSL --retry 3 --retry-delay 5 -o "${DATA_DIR}/Rfam.cm.gz" "${BASE_URL}/Rfam.cm.gz"
  gunzip -f "${DATA_DIR}/Rfam.cm.gz"
else
  echo "Already present, skipping: Rfam.cm"
fi

if [[ ! -s "${DATA_DIR}/Rfam.clanin" ]]; then
  curl -fSL --retry 3 --retry-delay 5 -o "${DATA_DIR}/Rfam.clanin" "${BASE_URL}/Rfam.clanin"
else
  echo "Already present, skipping: Rfam.clanin"
fi

test -s "${DATA_DIR}/Rfam.cm" || die "download appears incomplete: ${DATA_DIR}/Rfam.cm not found"
test -s "${DATA_DIR}/Rfam.clanin" || die "download appears incomplete: ${DATA_DIR}/Rfam.clanin not found"

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

echo "Indexing Rfam.cm with cmpress using ${RUNTIME} (${CONTAINER})"

case "$RUNTIME" in
  apptainer|singularity)
    "$RUNTIME" exec --bind "${DATA_DIR}:${DATA_DIR}" "docker://${CONTAINER}" \
      cmpress -F "${DATA_DIR}/Rfam.cm"
    ;;
  docker)
    docker run --rm -v "${DATA_DIR}:${DATA_DIR}" "${CONTAINER}" \
      cmpress -F "${DATA_DIR}/Rfam.cm"
    ;;
  *)
    die "unsupported runtime: ${RUNTIME}"
    ;;
esac

index_present || die "cmpress appears to have failed: ${DATA_DIR}/Rfam.cm.i1m not found"

echo "Rfam data ready and indexed at ${DATA_DIR}"
