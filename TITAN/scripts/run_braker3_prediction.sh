#!/usr/bin/env bash
set -euo pipefail

process_name="${1:?process name is required}"
genome="${2:?genome FASTA is required}"
task_cpus="${3:?task cpus is required}"
task_memory="${4:?task memory is required}"
has_long_reads="${5:?has_long_reads true/false is required}"

exec > >(tee -a braker3_run.log) 2>&1

BRAKER_DIR="/BRAKER-3.0.8"
PROTHINT_BIN="/ProtHint-2.6.0/bin"
GENEMARK_DIR="/GeneMark-ETP"
AUGUSTUS_CONFIG_DIR="/Augustus/config"
TSEBRA_BIN="/TSEBRA/bin"

log() {
    printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"
}

require_path() {
    local path="$1"
    local description="$2"
    if [[ ! -e "$path" ]]; then
        log "ERROR: missing ${description}: ${path}"
        exit 1
    fi
}

collect_files() {
    local directory="$1"
    local pattern="$2"
    local output_list="$3"
    local description="$4"

    if [[ ! -d "$directory" ]]; then
        log "ERROR: missing ${description} directory: ${directory}"
        exit 1
    fi

    find "$directory" -type f -name "$pattern" -print0 | sort -z | xargs -0 -r -n1 printf '%s\n' > "$output_list"
    if [[ ! -s "$output_list" ]]; then
        log "ERROR: no ${description} found in ${directory} with pattern ${pattern}"
        exit 1
    fi
}

copy_required_output() {
    local source="$1"
    local destination="$2"
    if [[ ! -s "$source" ]]; then
        log "ERROR: required BRAKER3 output is missing or empty: ${source}"
        log "Top-level files available after BRAKER3 failure/debug:"
        find . -maxdepth 2 -type f | sort | sed 's#^./#  #'
        exit 1
    fi
    cp "$source" "$destination"
}

log "Running BRAKER3/AUGUSTUS-GeneMark prediction"
log "Long-read BAM evidence enabled: ${has_long_reads}"
log "Work directory: ${PWD}"
log "Task resources: cpus=${task_cpus}; memory=${task_memory}"

require_path "${BRAKER_DIR}/scripts/braker.pl" "BRAKER braker.pl script"
require_path "${PROTHINT_BIN}" "ProtHint bin directory"
require_path "${GENEMARK_DIR}" "GeneMark-ETP directory"
require_path "${AUGUSTUS_CONFIG_DIR}" "AUGUSTUS config directory"
require_path "${TSEBRA_BIN}" "TSEBRA bin directory"
require_path "${genome}" "genome FASTA"

collect_files "bam_short" "*.bam" "bam_short.list" "short-read BAM evidence"
collect_files "protein_fastas" "*.fa" "protein_fastas.list" "cleaned protein FASTA evidence"

if [[ "${has_long_reads}" == "true" ]]; then
    collect_files "bam_long" "*.bam" "bam_long.list" "long-read BAM evidence"
    cat bam_short.list bam_long.list > bam_all.list
else
    : > bam_long.list
    cp bam_short.list bam_all.list
fi

bam_short_count=$(wc -l < bam_short.list | tr -d ' ')
bam_long_count=$(wc -l < bam_long.list | tr -d ' ')
bam_all_count=$(wc -l < bam_all.list | tr -d ' ')
protein_count=$(wc -l < protein_fastas.list | tr -d ' ')
bam=$(paste -sd, bam_all.list)
cleaned_proteins=$(paste -sd, protein_fastas.list)

{
    printf 'type\tpath\tsize_bytes\n'
    printf 'genome\t%s\t%s\n' "${genome}" "$(stat -c '%s' "${genome}")"
    while IFS= read -r file; do
        printf 'bam_short\t%s\t%s\n' "$file" "$(stat -c '%s' "$file")"
    done < bam_short.list
    while IFS= read -r file; do
        printf 'bam_long\t%s\t%s\n' "$file" "$(stat -c '%s' "$file")"
    done < bam_long.list
    while IFS= read -r file; do
        printf 'protein_fasta\t%s\t%s\n' "$file" "$(stat -c '%s' "$file")"
    done < protein_fastas.list
} > braker3_inputs.tsv

log "Short-read BAM files: ${bam_short_count}"
log "Long-read BAM files: ${bam_long_count}"
log "All BAM files: ${bam_all_count}"
log "Cleaned protein FASTA files: ${protein_count}"

braker_cmd=(
    "${BRAKER_DIR}/scripts/braker.pl"
    "--genome=${genome}"
    "--bam=${bam}"
    "--prot_seq=${cleaned_proteins}"
    "--threads=${task_cpus}"
    "--workingdir=${PWD}"
    "--softmasking"
    "--gff3"
    "--PROTHINT_PATH=${PROTHINT_BIN}/"
    "--GENEMARK_PATH=${GENEMARK_DIR}"
    "--AUGUSTUS_CONFIG_PATH=${AUGUSTUS_CONFIG_DIR}"
    "--TSEBRA_PATH=${TSEBRA_BIN}"
)

printf '%q ' "${braker_cmd[@]}" > braker3_command.txt
printf '\n' >> braker3_command.txt
log "Executing BRAKER3 command recorded in braker3_command.txt"

"${braker_cmd[@]}"

copy_required_output "Augustus/augustus.hints.gff3" "augustus.hints.gff3"
copy_required_output "GeneMark-ETP/genemark.gtf" "genemark.gtf"
copy_required_output "GeneMark-ETP/genemark_supported.gtf" "genemark_supported.gtf"
copy_required_output "braker.gff3" "braker.gff3"

{
    printf '"%s":\n' "${process_name}"
    printf '  braker: "3.0.8"\n'
    printf '  prothint: "2.6.0"\n'
    printf '  genemark_etp_path: "%s"\n' "${GENEMARK_DIR}"
    printf '  augustus_config_path: "%s"\n' "${AUGUSTUS_CONFIG_DIR}"
    printf '  tsebra_path: "%s"\n' "${TSEBRA_BIN}"
    printf '  bam_short_count: "%s"\n' "${bam_short_count}"
    printf '  bam_long_count: "%s"\n' "${bam_long_count}"
    printf '  bam_all_count: "%s"\n' "${bam_all_count}"
    printf '  protein_count: "%s"\n' "${protein_count}"
} > versions.yml

log "BRAKER3 completed successfully"
