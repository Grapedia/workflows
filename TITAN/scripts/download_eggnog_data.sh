#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
  cat <<'EOF'
Usage:
  scripts/download_eggnog_data.sh --data-dir PATH [--base-url URL]

Downloads the eggNOG-mapper reference database (eggnog.db, eggnog_proteins.dmnd
and eggnog.taxa.db) into --data-dir with plain curl/gunzip/tar. No eggNOG-mapper
installation or container runtime is required for this step.

The bundled download_eggnog_data.py script inside the eggnog-mapper 2.1.15
BioContainers image points at the retired host eggnogdb.embl.de, which no
longer resolves; this script fetches the same emapperdb-5.0.2 data version
straight from the current host instead.

Pass the resulting directory as --eggnog_data_dir to TITAN, together with
--run_eggnog_mapper true.

Options:
  --data-dir PATH  Target directory for the eggNOG-mapper database. Created if missing.
  --base-url URL   Base URL serving eggnog.db.gz, eggnog.taxa.tar.gz and
                   eggnog_proteins.dmnd.gz. Default: http://eggnog5.embl.de/download/emapperdb-5.0.2
  -h, --help       Show this help
EOF
}

die() {
  echo "ERROR: $*" >&2
  exit 1
}

DATA_DIR=""
BASE_URL="http://eggnog5.embl.de/download/emapperdb-5.0.2"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --data-dir) DATA_DIR="${2:-}"; shift 2 ;;
    --base-url) BASE_URL="${2:-}"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown option: $1. Use --help." ;;
  esac
done

[[ -n "$DATA_DIR" ]] || die "--data-dir is required"

command -v curl >/dev/null 2>&1 || die "curl is not available in PATH"
command -v gunzip >/dev/null 2>&1 || die "gunzip is not available in PATH"
command -v tar >/dev/null 2>&1 || die "tar is not available in PATH"

mkdir -p "$DATA_DIR"
DATA_DIR="$(cd "$DATA_DIR" && pwd -P)"

echo "Downloading eggNOG-mapper database into ${DATA_DIR} from ${BASE_URL}"

fetch_and_unpack() {
  local remote_name="$1"
  local extract_mode="$2"
  local target_name="$3"

  if [[ -s "${DATA_DIR}/${target_name}" ]]; then
    echo "Already present, skipping: ${target_name}"
    return 0
  fi

  echo "Fetching ${remote_name}..."
  curl -fSL --retry 3 --retry-delay 5 -o "${DATA_DIR}/${remote_name}" "${BASE_URL}/${remote_name}"

  case "$extract_mode" in
    gunzip)
      gunzip -f "${DATA_DIR}/${remote_name}"
      ;;
    tar)
      tar -xzf "${DATA_DIR}/${remote_name}" -C "${DATA_DIR}"
      rm -f "${DATA_DIR}/${remote_name}"
      ;;
  esac
}

fetch_and_unpack "eggnog.db.gz" gunzip "eggnog.db"
fetch_and_unpack "eggnog.taxa.tar.gz" tar "eggnog.taxa.db"
fetch_and_unpack "eggnog_proteins.dmnd.gz" gunzip "eggnog_proteins.dmnd"

test -s "${DATA_DIR}/eggnog.db" || die "download appears incomplete: ${DATA_DIR}/eggnog.db not found"
test -s "${DATA_DIR}/eggnog_proteins.dmnd" || die "download appears incomplete: ${DATA_DIR}/eggnog_proteins.dmnd not found"
test -s "${DATA_DIR}/eggnog.taxa.db" || die "download appears incomplete: ${DATA_DIR}/eggnog.taxa.db not found"

echo "eggNOG-mapper database ready at ${DATA_DIR}"
