#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
  cat <<'EOF'
Usage:
  scripts/download_interproscan_data.sh --data-dir PATH [--version VERSION] [--base-url URL]

Downloads the InterProScan member database data-only bundle (pfam, cdd,
gene3d, panther, etc.) into --data-dir with plain curl/tar, verified against
the upstream .md5 checksum. No InterProScan installation or container
runtime is required for this step.

The resulting directory is meant to be bind-mounted to /opt/interproscan/data
inside the interpro/interproscan container, matching the official Docker
usage pattern (`-v .../interproscan-<version>/data:/opt/interproscan/data`).

Pass the resulting directory as --interproscan_data_dir to TITAN, together
with --run_interproscan true.

Options:
  --data-dir PATH  Target directory for the InterProScan member database data.
                   Created if missing. Populated directly with the pfam/,
                   cdd/, gene3d/, ... subdirectories (no extra nesting).
  --version VER    InterProScan version. Must match container_interproscan.
                   Default: 5.78-109.0
  --base-url URL   Base URL serving the versioned interproscan-data-<version>.tar.gz
                   and its .md5 file.
                   Default: https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5
  -h, --help       Show this help
EOF
}

die() {
  echo "ERROR: $*" >&2
  exit 1
}

DATA_DIR=""
VERSION="5.78-109.0"
BASE_URL="https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --data-dir) DATA_DIR="${2:-}"; shift 2 ;;
    --version) VERSION="${2:-}"; shift 2 ;;
    --base-url) BASE_URL="${2:-}"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown option: $1. Use --help." ;;
  esac
done

[[ -n "$DATA_DIR" ]] || die "--data-dir is required"

command -v curl >/dev/null 2>&1 || die "curl is not available in PATH"
command -v tar >/dev/null 2>&1 || die "tar is not available in PATH"
command -v md5sum >/dev/null 2>&1 || die "md5sum is not available in PATH"

mkdir -p "$DATA_DIR"
DATA_DIR="$(cd "$DATA_DIR" && pwd -P)"

if [[ -d "${DATA_DIR}/pfam" ]] && find "${DATA_DIR}/pfam" -type f -print -quit | grep -q .; then
  echo "InterProScan member database data already present, skipping: ${DATA_DIR}/pfam"
  exit 0
fi

ARCHIVE="interproscan-data-${VERSION}.tar.gz"
ARCHIVE_URL="${BASE_URL}/${VERSION}/alt/${ARCHIVE}"
MD5_URL="${ARCHIVE_URL}.md5"

echo "Downloading InterProScan ${VERSION} member database data into ${DATA_DIR} from ${ARCHIVE_URL}"
echo "This is a large download (several GB compressed); it may take a while."

WORK_DIR="$(mktemp -d)"
trap 'rm -rf "$WORK_DIR"' EXIT

curl -fSL --retry 3 --retry-delay 5 -o "${WORK_DIR}/${ARCHIVE}" "${ARCHIVE_URL}"
curl -fSL --retry 3 --retry-delay 5 -o "${WORK_DIR}/${ARCHIVE}.md5" "${MD5_URL}"

(
  cd "$WORK_DIR"
  md5sum -c "${ARCHIVE}.md5"
) || die "checksum verification failed for ${ARCHIVE}"

# The archive's top-level entries are interproscan-<version>/data/<member-db>/...;
# strip both components so $DATA_DIR ends up holding pfam/, cdd/, gene3d/, ...
# directly, matching the bind-mount target /opt/interproscan/data.
tar -xzf "${WORK_DIR}/${ARCHIVE}" --strip-components=2 -C "$DATA_DIR"

test -d "${DATA_DIR}/pfam" || die "download appears incomplete: ${DATA_DIR}/pfam not found"

echo "InterProScan member database data ready at ${DATA_DIR}"
