#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 20 ]]; then
  echo "Usage: run_aegis_merge.sh <mode> <masked_genome> <aegis_version> <aegis_container> <output_prefix> <liftoff_gff3> <augustus_gff3> <genemark_gtf> <egapx_gff3> <star_stringtie_stranded_default_gtf> <star_stringtie_stranded_alt_gtf> <star_psiclass_stranded_gtf> <long_reads_default_gtf> <long_reads_alt_gtf> <flair_isoforms_gtf> <star_psiclass_unstranded_gtf> <star_stringtie_unstranded_default_gtf> <star_stringtie_unstranded_alt_gtf> <helixer_gff3> <gene_id_prefix>" >&2
  exit 2
fi

mode="$1"
masked_genome="$2"
aegis_version="$3"
aegis_container="$4"
output_prefix="$5"
liftoff_gff3="$6"
augustus_gff3="$7"
genemark_gtf="$8"
egapx_gff3="$9"
star_stringtie_stranded_default_gtf="${10}"
star_stringtie_stranded_alt_gtf="${11}"
star_psiclass_stranded_gtf="${12}"
long_reads_default_gtf="${13}"
long_reads_alt_gtf="${14}"
flair_isoforms_gtf="${15}"
star_psiclass_unstranded_gtf="${16}"
star_stringtie_unstranded_default_gtf="${17}"
star_stringtie_unstranded_alt_gtf="${18}"
helixer_gff3="${19}"
gene_id_prefix="${20}"

if [[ "$mode" != "short_reads" && "$mode" != "short_and_long_reads" ]]; then
  echo "AEGIS mode must be 'short_reads' or 'short_and_long_reads', got: $mode" >&2
  exit 1
fi

date_stamp="$(date "+%Y-%m-%d %H:%M:%S")"
echo "[$date_stamp] Running AEGIS merge on ${mode} evidence"

has_feature_records() {
  local file="$1"
  awk 'NF && $0 !~ /^#/ { found = 1; exit } END { exit found ? 0 : 1 }' "$file"
}

require_nonempty_evidence() {
  local name="$1"
  local file="$2"
  if [[ ! -f "$file" ]]; then
    echo "Required AEGIS evidence is missing (${name}): ${file}" >&2
    exit 1
  fi
  if [[ ! -s "$file" || ! -r "$file" ]]; then
    echo "Required AEGIS evidence is empty or unreadable (${name}): ${file}" >&2
    exit 1
  fi
  if ! has_feature_records "$file"; then
    echo "Required AEGIS evidence has no feature records (${name}): ${file}" >&2
    exit 1
  fi
}

include_optional_evidence() {
  local name="$1"
  local file="$2"
  if [[ ! -f "$file" || ! -s "$file" ]]; then
    printf '%s\t%s\tfalse\tfalse\t0\n' "$name" "$file" >> aegis_inputs.tsv
    return 0
  fi
  if has_feature_records "$file"; then
    merge_inputs+=("$file")
    printf '%s\t%s\tfalse\ttrue\t%s\n' "$name" "$file" "$(stat -c '%s' "$file")" >> aegis_inputs.tsv
  else
    printf '%s\t%s\tfalse\tfalse\t%s\n' "$name" "$file" "$(stat -c '%s' "$file")" >> aegis_inputs.tsv
  fi
}

include_required_evidence() {
  local name="$1"
  local file="$2"
  require_nonempty_evidence "$name" "$file"
  merge_inputs+=("$file")
  printf '%s\t%s\ttrue\ttrue\t%s\n' "$name" "$file" "$(stat -c '%s' "$file")" >> aegis_inputs.tsv
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

merge_inputs=()
printf 'name\tpath\trequired\tincluded\tsize_bytes\n' > aegis_inputs.tsv

include_required_evidence "liftoff_gff3" "$liftoff_gff3"
include_required_evidence "augustus_gff3" "$augustus_gff3"
include_required_evidence "genemark_gtf" "$genemark_gtf"
include_required_evidence "egapx_gff3" "$egapx_gff3"
include_required_evidence "star_stringtie_stranded_default_gtf" "$star_stringtie_stranded_default_gtf"
include_required_evidence "star_stringtie_stranded_alt_gtf" "$star_stringtie_stranded_alt_gtf"
include_required_evidence "star_psiclass_stranded_gtf" "$star_psiclass_stranded_gtf"

if [[ "$mode" == "short_and_long_reads" ]]; then
  include_required_evidence "long_reads_default_gtf" "$long_reads_default_gtf"
  include_required_evidence "long_reads_alt_gtf" "$long_reads_alt_gtf"
else
  include_optional_evidence "long_reads_default_gtf" "$long_reads_default_gtf"
  include_optional_evidence "long_reads_alt_gtf" "$long_reads_alt_gtf"
fi

include_optional_evidence "star_psiclass_unstranded_gtf" "$star_psiclass_unstranded_gtf"
include_optional_evidence "star_stringtie_unstranded_default_gtf" "$star_stringtie_unstranded_default_gtf"
include_optional_evidence "star_stringtie_unstranded_alt_gtf" "$star_stringtie_unstranded_alt_gtf"
include_optional_evidence "flair_isoforms_gtf" "$flair_isoforms_gtf"
include_optional_evidence "helixer_gff3" "$helixer_gff3"

if [[ "${#merge_inputs[@]}" -eq 0 ]]; then
  echo "No non-empty AEGIS evidence files were provided" >&2
  exit 1
fi

printf "[%s] AEGIS merge inputs:\n" "$date_stamp"
printf '  %s\n' "${merge_inputs[@]}"

aegis_cmd=(/opt/conda/envs/bio_env/bin/python -m aegis)

echo "[$date_stamp] AEGIS merge"
"${aegis_cmd[@]}" merge -d aegis_merge -o "$output_prefix" "${merge_inputs[@]}"
if [[ ! -s "aegis_merge/${output_prefix}.gff3" ]]; then
  echo "AEGIS merge did not produce ${output_prefix}.gff3" >&2
  exit 1
fi

echo "[$date_stamp] AEGIS rename (prefix: ${gene_id_prefix})"
"${aegis_cmd[@]}" rename \
  -a "$output_prefix" \
  -d aegis_rename \
  --prefix "$gene_id_prefix" \
  --gene-id-correspondences \
  "aegis_merge/${output_prefix}.gff3"
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

printf '"%s":\n  aegis: "%s"\n  aegis_container: "%s"\n  mode: "%s"\n  evidence_count: "%s"\n  gene_id_prefix: "%s"\n' \
  "${NXF_TASK_PROCESS:-aegis}" "$aegis_version" "$aegis_container" "$mode" "${#merge_inputs[@]}" "$gene_id_prefix" > versions.yml
