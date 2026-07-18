#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"
PYTHON_BIN="${PYTHON_BIN:-$(command -v python3.11 || command -v python3)}"

run_step() {
  printf '\n==> %s\n' "$*"
  "$@"
}

expect_failure() {
  local name="$1"
  local expected="$2"
  shift 2
  local output_file
  output_file="$(mktemp)"

  printf '\n==> %s\n' "$name"
  if "$@" >"$output_file" 2>&1; then
    cat "$output_file"
    rm -f "$output_file"
    printf 'ERROR: %s: expected failure but command succeeded\n' "$name" >&2
    return 1
  fi

  if ! grep -Fq "$expected" "$output_file"; then
    cat "$output_file"
    rm -f "$output_file"
    printf 'ERROR: %s: expected message fragment not found: %s\n' "$name" "$expected" >&2
    return 1
  fi

  rm -f "$output_file"
}

run_step "$PYTHON_BIN" scripts/validate_container_pins.py
run_step "$PYTHON_BIN" scripts/validate_nextflow_quality.py
run_step "$PYTHON_BIN" scripts/validate_profiles.py
run_step "$PYTHON_BIN" scripts/validate_minimal_test_data.py
run_step "$PYTHON_BIN" scripts/test_validate_inputs.py
run_step "$PYTHON_BIN" scripts/test_validate_final_annotation.py
run_step "$PYTHON_BIN" scripts/test_download_sra_fastq.py
run_step "$PYTHON_BIN" scripts/test_shared_shell_scripts.py
run_step "$PYTHON_BIN" scripts/test_clean_liftoff_gff3_for_agat.py
run_step "$PYTHON_BIN" scripts/test_trnascan_to_gff3.py
run_step "$PYTHON_BIN" scripts/test_rfam_tblout_to_gff3.py
run_step "$PYTHON_BIN" scripts/test_build_lncrna_candidates.py
run_step "$PYTHON_BIN" scripts/test_make_mikado_list.py
run_step "$PYTHON_BIN" scripts/test_compare_final_annotations.py
run_step bash scripts/test_download_cpat_plant_lncpipe.sh

run_step nextflow config -profile test
run_step nextflow run main.nf -profile test -stub-run -ansi-log false
run_step "$PYTHON_BIN" -c 'import json; report=json.load(open("test-results/validation/final_annotation_validation.json")); assert report["status"] == "pass", report'
run_step nextflow run main.nf -profile test -stub-run -ansi-log false \
  --RNAseq_samplesheet test-data/minimal/valid/rnaseq_samplesheet_short_reads_only.csv \
  --output_dir test-results/short-reads-only
run_step "$PYTHON_BIN" -c 'import json; manifest=json.load(open("test-results/short-reads-only/provenance/evidence_manifest.json")); assert manifest["has_long_reads"] is False, manifest'
run_step nextflow run main.nf -profile test -stub-run -ansi-log false \
  --RNAseq_samplesheet test-data/minimal/valid/rnaseq_samplesheet_single_sample.csv \
  --output_dir test-results/single-sample
run_step "$PYTHON_BIN" -c 'import json; report=json.load(open("test-results/single-sample/validation/final_annotation_validation.json")); assert report["status"] == "pass", report'

expect_failure \
  "Nextflow rejects invalid RNA-seq samplesheet before heavy execution" \
  "library_layout must be one of" \
  nextflow run main.nf -profile test \
    --RNAseq_samplesheet test-data/minimal/invalid/rnaseq_bad_layout.csv \
    -stub-run -ansi-log false

printf '\nTITAN quick tests OK\n'
