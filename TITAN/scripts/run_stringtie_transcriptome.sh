#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 7 ]]; then
  echo "Usage: run_stringtie_transcriptome.sh <stringtie_script> <alt_script> <threads> <bam> <read_type> <default_gtf> <alt_gtf>" >&2
  exit 2
fi

stringtie_script="$1"
alt_script="$2"
threads="$3"
bam_file="$4"
read_type="$5"
default_gtf="$6"
alt_gtf="$7"

date_stamp="$(date "+%Y-%m-%d %H:%M:%S")"

cmd=(bash "$stringtie_script" -t "$threads" -o "$default_gtf" -b "$bam_file" -r "$read_type")
echo "[$date_stamp] Executing: ${cmd[*]}"
"${cmd[@]}"

cmd=(bash "$alt_script" -t "$threads" -o "$alt_gtf" -b "$bam_file" -r "$read_type")
echo "[$date_stamp] Executing: ${cmd[*]}"
"${cmd[@]}"
