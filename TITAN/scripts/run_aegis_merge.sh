#!/usr/bin/env bash
set -euo pipefail

# NOTE: This script used to fan out to `aegis merge` across ~13 raw evidence
# GFF/GTF sources (liftoff, BRAKER, egapx, helixer, and every transcript
# assembly). That step is now skipped: `mikado_pick` already consolidates
# those same sources into a single non-redundant, per-locus-scored gene set
# (numeric source_score/reference weighting in mikado_configuration.yaml,
# built from the same evidence), which is a more faithful reconciliation
# than aegis merge's whole-gene overlap-threshold exclusion. AEGIS is now
# only used here to rename IDs (Vitvi prefix) and tidy the Mikado output.

if [[ $# -ne 6 ]]; then
  echo "Usage: run_aegis_merge.sh <masked_genome> <aegis_version> <aegis_container> <output_prefix> <mikado_gff3> <gene_id_prefix>" >&2
  exit 2
fi

masked_genome="$1"
aegis_version="$2"
aegis_container="$3"
output_prefix="$4"
mikado_gff3="$5"
gene_id_prefix="$6"

date_stamp="$(date "+%Y-%m-%d %H:%M:%S")"
echo "[$date_stamp] Finalizing annotation from mikado_pick output"

has_feature_records() {
  local file="$1"
  awk 'NF && $0 !~ /^#/ { found = 1; exit } END { exit found ? 0 : 1 }' "$file"
}

require_nonempty_file() {
  local name="$1"
  local file="$2"
  if [[ ! -f "$file" || ! -s "$file" ]]; then
    echo "Required AEGIS input is missing or empty (${name}): ${file}" >&2
    exit 1
  fi
}

require_nonempty_file "masked genome" "$masked_genome"
require_nonempty_file "mikado_pick annotation" "$mikado_gff3"
if ! has_feature_records "$mikado_gff3"; then
  echo "mikado_pick annotation has no feature records (${mikado_gff3}); check --run_mikado and upstream evidence" >&2
  exit 1
fi

printf 'name\tpath\trequired\tincluded\tsize_bytes\n' > aegis_inputs.tsv
printf 'mikado_pick\t%s\ttrue\ttrue\t%s\n' "$mikado_gff3" "$(stat -c '%s' "$mikado_gff3")" >> aegis_inputs.tsv

aegis_cmd=(/opt/conda/envs/bio_env/bin/python -m aegis)

echo "[$date_stamp] AEGIS rename (prefix: ${gene_id_prefix})"
"${aegis_cmd[@]}" rename \
  -a "$output_prefix" \
  -d aegis_rename \
  --prefix "$gene_id_prefix" \
  --gene-id-correspondences \
  "$mikado_gff3"
if [[ ! -s "aegis_rename/${output_prefix}_renamed.gff3" ]]; then
  echo "AEGIS rename did not produce ${output_prefix}_renamed.gff3" >&2
  exit 1
fi

echo "[$date_stamp] AEGIS tidy"
"${aegis_cmd[@]}" tidy \
  -a "${output_prefix}_renamed" \
  -d aegis_tidy \
  --standard-features \
  "aegis_rename/${output_prefix}_renamed.gff3"
if [[ ! -s "aegis_tidy/${output_prefix}_renamed_tidy.gff3" ]]; then
  echo "AEGIS tidy did not produce ${output_prefix}_renamed_tidy.gff3" >&2
  exit 1
fi
cp "aegis_tidy/${output_prefix}_renamed_tidy.gff3" final_annotation.gff3

echo "[$date_stamp] AEGIS extract"
"${aegis_cmd[@]}" extract -f protein -m all -d aegis_proteins_all final_annotation.gff3 "$masked_genome"
"${aegis_cmd[@]}" extract -f protein -m main -d aegis_proteins_main final_annotation.gff3 "$masked_genome"

all_proteins="$(find aegis_proteins_all -type f -name '*proteins*all*.fasta' -print -quit 2>/dev/null || true)"
main_proteins="$(find aegis_proteins_main -type f -name '*proteins*main*.fasta' -print -quit 2>/dev/null || true)"

if [[ -z "$all_proteins" ]]; then
  echo "AEGIS did not produce an all-proteins FASTA" >&2
  exit 1
fi

if [[ -z "$main_proteins" ]]; then
  echo "AEGIS did not produce a main-proteins FASTA" >&2
  exit 1
fi

cp "$all_proteins" final_annotation_proteins_all.fasta
cp "$main_proteins" final_annotation_proteins_main.fasta
test -s final_annotation.gff3
test -s final_annotation_proteins_all.fasta
test -s final_annotation_proteins_main.fasta

printf '"%s":\n  aegis: "%s"\n  aegis_container: "%s"\n  source: "mikado_pick"\n  gene_id_prefix: "%s"\n' \
  "${NXF_TASK_PROCESS:-aegis}" "$aegis_version" "$aegis_container" "$gene_id_prefix" > versions.yml
