#!/usr/bin/env bash
set -Eeuo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_DIR"' EXIT
cd "$ROOT_DIR"

"${ROOT_DIR}/scripts/download_cpat_plant_lncpipe.sh" \
  --model-dir "${TMP_DIR}/model" \
  --base-url "file://${ROOT_DIR}/resources/cpat_plant_lncpipe"

test -s "${TMP_DIR}/model/Plant_Hexamer.tsv"
test -s "${TMP_DIR}/model/Plant.logit.RData"
sha256sum -c "resources/cpat_plant_lncpipe/SHA256SUMS"

second_run_output="$("${ROOT_DIR}/scripts/download_cpat_plant_lncpipe.sh" \
  --model-dir "${TMP_DIR}/model" \
  --base-url "file://${ROOT_DIR}/resources/cpat_plant_lncpipe")"
grep -q "Already present" <<<"$second_run_output"

echo "download_cpat_plant_lncpipe tests OK"
