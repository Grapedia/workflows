#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
  cat <<'EOF'
Usage:
  scripts/download_cpat_plant_lncpipe.sh --model-dir PATH [--base-url URL]

Downloads the Plant-LncPipe CPAT-plant model files into --model-dir:
Plant_Hexamer.tsv and Plant.logit.RData. Existing valid files are left in
place. Files are verified with pinned SHA-256 checksums.

Options:
  --model-dir PATH  Target directory for the CPAT-plant model files.
  --base-url URL    Base raw GitHub URL containing Plant_Hexamer.tsv and
                    Plant.logit.RData.
                    Default: https://raw.githubusercontent.com/xuechantian/Plant-LncRNA-pipline/master/Model
  -h, --help        Show this help
EOF
}

die() {
  echo "ERROR: $*" >&2
  exit 1
}

MODEL_DIR=""
BASE_URL="https://raw.githubusercontent.com/xuechantian/Plant-LncRNA-pipline/master/Model"
HEXAMER_SHA256="db99bc9c5f1619326019e64dcabfa89be7c4dac881357b81b9ebf265295bb60c"
LOGIT_SHA256="774b612fd0a017a282367d59471838bf6328d99215f1fd9ad1cbf63457e1071e"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --model-dir) MODEL_DIR="${2:-}"; shift 2 ;;
    --base-url) BASE_URL="${2:-}"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown option: $1. Use --help." ;;
  esac
done

[[ -n "$MODEL_DIR" ]] || die "--model-dir is required"

if ! command -v curl >/dev/null 2>&1 && ! command -v python3 >/dev/null 2>&1; then
  die "neither curl nor python3 is available in PATH"
fi
command -v sha256sum >/dev/null 2>&1 || die "sha256sum is not available in PATH"

mkdir -p "$MODEL_DIR"
MODEL_DIR="$(cd "$MODEL_DIR" && pwd -P)"

verify_file() {
  local path="$1"
  local expected="$2"
  [[ -s "$path" ]] || return 1
  printf '%s  %s\n' "$expected" "$path" | sha256sum -c - >/dev/null 2>&1
}

fetch_file() {
  local filename="$1"
  local expected="$2"
  local destination="${MODEL_DIR}/${filename}"

  if verify_file "$destination" "$expected"; then
    echo "Already present and checksum-valid, skipping: ${destination}"
    return 0
  fi

  echo "Fetching ${filename} from ${BASE_URL}"
  local tmp
  tmp="$(mktemp "${MODEL_DIR}/.${filename}.XXXXXX")"
  if command -v curl >/dev/null 2>&1; then
    curl -fSL --retry 3 --retry-delay 5 -o "$tmp" "${BASE_URL}/${filename}"
  else
    python3 - "$tmp" "${BASE_URL}/${filename}" <<'PY'
import shutil
import sys
from urllib.request import urlopen

destination, url = sys.argv[1:3]
with urlopen(url, timeout=60) as response, open(destination, "wb") as handle:
    shutil.copyfileobj(response, handle)
PY
  fi
  printf '%s  %s\n' "$expected" "$tmp" | sha256sum -c - >/dev/null
  mv -f "$tmp" "$destination"
}

fetch_file "Plant_Hexamer.tsv" "$HEXAMER_SHA256"
fetch_file "Plant.logit.RData" "$LOGIT_SHA256"

cat > "${MODEL_DIR}/SHA256SUMS" <<EOF
${HEXAMER_SHA256}  ${MODEL_DIR}/Plant_Hexamer.tsv
${LOGIT_SHA256}  ${MODEL_DIR}/Plant.logit.RData
EOF

echo "CPAT-plant Plant-LncPipe model ready at ${MODEL_DIR}"
