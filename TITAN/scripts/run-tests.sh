#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

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

run_step python3 scripts/validate_container_pins.py
run_step python3 scripts/validate_nextflow_quality.py
run_step python3 scripts/validate_profiles.py
run_step python3 scripts/validate_minimal_test_data.py
run_step python3 scripts/test_validate_inputs.py
run_step python3 scripts/test_validate_final_annotation.py
run_step python3 scripts/test_download_sra_fastq.py
run_step python3 scripts/test_shared_shell_scripts.py

run_step nextflow config -profile test
run_step nextflow run main.nf -profile test -stub-run -ansi-log false
run_step python3 -c 'import json; report=json.load(open("test-results/validation/final_annotation_validation.json")); assert report["status"] == "pass", report'
run_step nextflow run main.nf -profile test -stub-run -ansi-log false \
  --RNAseq_samplesheet test-data/minimal/valid/rnaseq_samplesheet_short_reads_only.csv \
  --output_dir test-results/short-reads-only
run_step python3 -c 'import json; manifest=json.load(open("test-results/short-reads-only/provenance/evidence_manifest.json")); assert manifest["has_long_reads"] is False, manifest'
run_step nextflow run main.nf -profile test -stub-run -ansi-log false \
  --RNAseq_samplesheet test-data/minimal/valid/rnaseq_samplesheet_single_sample.csv \
  --output_dir test-results/single-sample
run_step python3 -c 'import json; report=json.load(open("test-results/single-sample/validation/final_annotation_validation.json")); assert report["status"] == "pass", report'

expect_failure \
  "Nextflow rejects invalid RNA-seq samplesheet before heavy execution" \
  "library_layout must be one of" \
  nextflow run main.nf -profile test \
    --RNAseq_samplesheet test-data/minimal/invalid/rnaseq_bad_layout.csv \
    -stub-run -ansi-log false

printf '\nTITAN quick tests OK\n'
