#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 6 ]]; then
  echo "Usage: run_aegis_merge.sh <mode> <masked_genome> <aegis_version> <aegis_container> <output_prefix> <evidence...>" >&2
  exit 2
fi

mode="$1"
masked_genome="$2"
aegis_version="$3"
aegis_container="$4"
output_prefix="$5"
shift 5

date_stamp="$(date "+%Y-%m-%d %H:%M:%S")"
echo "[$date_stamp] Running AEGIS merge on ${mode} evidence"

merge_inputs=()
for evidence in "$@"; do
  if [[ -s "$evidence" ]]; then
    merge_inputs+=("$evidence")
  fi
done

if [[ "${#merge_inputs[@]}" -eq 0 ]]; then
  echo "No non-empty AEGIS evidence files were provided" >&2
  exit 1
fi

printf "[%s] AEGIS merge inputs:\n" "$date_stamp"
printf '  %s\n' "${merge_inputs[@]}"

aegis merge -d aegis_merge -o "$output_prefix" "${merge_inputs[@]}"
cp "aegis_merge/${output_prefix}.gff3" final_annotation.gff3

aegis extract -f protein -m all -d aegis_proteins_all final_annotation.gff3 "$masked_genome"
aegis extract -f protein -m main -d aegis_proteins_main final_annotation.gff3 "$masked_genome"

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

printf '"%s":\n  aegis: "%s"\n  aegis_container: "%s"\n' \
  "${NXF_TASK_PROCESS:-aegis}" "$aegis_version" "$aegis_container" > versions.yml
