#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
  cat <<'EOF'
Usage:
  scripts/download_omark_data.sh --data-dir PATH [--base-url URL]

Downloads the OMAmer reference database (omamer.h5) into --data-dir with
plain curl. No OMArk/OMAmer installation or container runtime is required for
this step.

Defaults to the full LUCA.h5 database (whole OMA database), as recommended by
the OMArk README over a taxonomically-restricted database: a narrower
database (e.g. Viridiplantae-only) limits OMArk's ability to detect
contamination or sequences from outside that range, which is one of the two
things OMArk is used for in TITAN (the other being completeness, alongside
BUSCO).

Pass the resulting directory as --omark_data_dir to TITAN, together with
--run_omark true.

Options:
  --data-dir PATH  Target directory for omamer.h5. Created if missing.
  --base-url URL   URL of the OMAmer .h5 database file to download.
                   Default: https://omabrowser.org/All/LUCA.h5
  -h, --help       Show this help
EOF
}

die() {
  echo "ERROR: $*" >&2
  exit 1
}

DATA_DIR=""
BASE_URL="https://omabrowser.org/All/LUCA.h5"

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

mkdir -p "$DATA_DIR"
DATA_DIR="$(cd "$DATA_DIR" && pwd -P)"

if [[ -s "${DATA_DIR}/omamer.h5" ]]; then
  echo "OMAmer database already present, skipping: ${DATA_DIR}/omamer.h5"
  exit 0
fi

echo "Downloading OMAmer database into ${DATA_DIR} from ${BASE_URL}"
echo "This is a large download (multi-GB); it may take a while."

TMP_FILE="$(mktemp "${DATA_DIR}/.omamer.h5.XXXXXX")"
curl -fSL --retry 3 --retry-delay 5 -o "$TMP_FILE" "$BASE_URL"
mv -f "$TMP_FILE" "${DATA_DIR}/omamer.h5"

test -s "${DATA_DIR}/omamer.h5" || die "download appears incomplete: ${DATA_DIR}/omamer.h5 not found"

echo "OMAmer database ready at ${DATA_DIR}/omamer.h5"
