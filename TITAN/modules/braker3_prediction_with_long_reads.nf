// AUGUSTUS can be run directly using the pipeline BRAKER3
process braker3_prediction_with_long_reads {
  label 'process_prediction'

  tag "Executing BRAKER3/AUGUSTUS-Genemark prediction"
  container params.container_braker3
  publishDir "${params.output_dir}", mode: 'copy', saveAs: { filename ->
    if (filename in ['augustus.hints.gff3', 'genemark.gtf', 'genemark_supported.gtf', 'braker.gff3', 'versions.yml']) {
      return filename
    }
    return null
  }
  publishDir "${params.output_dir}/intermediate_files/braker3", mode: 'copy', enabled: params.publish_intermediates, saveAs: { filename ->
    if (filename in ['braker3_run.log', 'braker3_command.txt', 'braker3_inputs.tsv', 'braker.log', 'errors']) {
      return filename
    }
    return null
  }

  input:
    path(genome)
    path(protein_fastas, stageAs: "protein_fastas/*")
    path(bam_short, stageAs: "bam_short/*")
    path(bam_long, stageAs: "bam_long/*")
    path(clean_protein_script)

  output:
    path "augustus.hints.gff3", emit: augustus_gff
    path "genemark.gtf", emit: genemark_gtf
    path "genemark_supported.gtf", emit: genemark_supported_gtf
    path "braker.gff3", emit: braker_gff
    path "braker3_run.log", emit: debug_log
    path "braker3_command.txt", emit: command_log
    path "braker3_inputs.tsv", emit: input_manifest
    path "braker.log", optional: true, emit: braker_log
    path "errors", optional: true, emit: errors
    path "versions.yml", emit: versions



  script:
    """
    set -euo pipefail
    exec > >(tee -a braker3_run.log) 2>&1

    BRAKER_DIR="/BRAKER-3.0.8"
    PROTHINT_BIN="/ProtHint-2.6.0/bin"
    GENEMARK_DIR="/GeneMark-ETP"
    AUGUSTUS_CONFIG_DIR="/Augustus/config"
    TSEBRA_BIN="/TSEBRA/bin"

    log() {
        printf '[%s] %s\\n' "\$(date '+%Y-%m-%d %H:%M:%S')" "\$*"
    }

    require_path() {
        local path="\$1"
        local description="\$2"
        if [[ ! -e "\$path" ]]; then
            log "ERROR: missing \${description}: \${path}"
            exit 1
        fi
    }

    collect_files() {
        local directory="\$1"
        local pattern="\$2"
        local output_list="\$3"
        local description="\$4"

        if [[ ! -d "\$directory" ]]; then
            log "ERROR: missing \${description} directory: \${directory}"
            exit 1
        fi

        find "\$directory" -type f -name "\$pattern" -print0 | sort -z | xargs -0 -r -n1 printf '%s\\n' > "\$output_list"
        if [[ ! -s "\$output_list" ]]; then
            log "ERROR: no \${description} found in \${directory} with pattern \${pattern}"
            exit 1
        fi
    }

    copy_required_output() {
        local source="\$1"
        local destination="\$2"
        if [[ ! -s "\$source" ]]; then
            log "ERROR: required BRAKER3 output is missing or empty: \${source}"
            log "Top-level files available after BRAKER3 failure/debug:"
            find . -maxdepth 2 -type f | sort | sed 's#^./#  #'
            exit 1
        fi
        cp "\$source" "\$destination"
    }

    log "Running BRAKER3/AUGUSTUS-GeneMark prediction with short- and long-read BAM evidence"
    log "Work directory: \${PWD}"
    log "Task resources: cpus=${task.cpus}; memory=${task.memory ?: 'not_set'}"

    require_path "\${BRAKER_DIR}/scripts/braker.pl" "BRAKER braker.pl script"
    require_path "\${PROTHINT_BIN}" "ProtHint bin directory"
    require_path "\${GENEMARK_DIR}" "GeneMark-ETP directory"
    require_path "\${AUGUSTUS_CONFIG_DIR}" "AUGUSTUS config directory"
    require_path "\${TSEBRA_BIN}" "TSEBRA bin directory"
    require_path "${clean_protein_script}" "protein cleanup script"
    require_path "${genome}" "genome FASTA"

    collect_files "bam_short" "*.bam" "bam_short.list" "short-read BAM evidence"
    collect_files "bam_long" "*.bam" "bam_long.list" "long-read BAM evidence"
    collect_files "protein_fastas" "*" "protein_fastas.list" "protein FASTA evidence"
    cat bam_short.list bam_long.list > bam_all.list

    bam_short_count=\$(wc -l < bam_short.list | tr -d ' ')
    bam_long_count=\$(wc -l < bam_long.list | tr -d ' ')
    bam_all_count=\$(wc -l < bam_all.list | tr -d ' ')
    protein_count=\$(wc -l < protein_fastas.list | tr -d ' ')
    bam=\$(paste -sd, bam_all.list)

    {
        printf 'type\\tpath\\tsize_bytes\\n'
        printf 'genome\\t%s\\t%s\\n' "${genome}" "\$(stat -c '%s' "${genome}")"
        while IFS= read -r file; do
            printf 'bam_short\\t%s\\t%s\\n' "\$file" "\$(stat -c '%s' "\$file")"
        done < bam_short.list
        while IFS= read -r file; do
            printf 'bam_long\\t%s\\t%s\\n' "\$file" "\$(stat -c '%s' "\$file")"
        done < bam_long.list
        while IFS= read -r file; do
            printf 'protein_fasta\\t%s\\t%s\\n' "\$file" "\$(stat -c '%s' "\$file")"
        done < protein_fastas.list
    } > braker3_inputs.tsv

    log "Short-read BAM files: \${bam_short_count}"
    log "Long-read BAM files: \${bam_long_count}"
    log "All BAM files: \${bam_all_count}"
    log "Protein FASTA files: \${protein_count}"
    log "Cleaning protein FASTA headers for BRAKER3"

    : > cleaned_proteins.list
    while IFS= read -r protein; do
        stem=\$(basename "\$protein")
        stem="\${stem%.*}"
        safe_stem=\$(printf '%s' "\$stem" | tr -c 'A-Za-z0-9_.-' '_')
        cleaned="\${safe_stem}.cleaned.fasta"
        python3 "${clean_protein_script}" "\$protein" "\$cleaned" "\$stem"
        if [[ ! -s "\$cleaned" ]]; then
            log "ERROR: cleaned protein FASTA is missing or empty: \${cleaned}"
            exit 1
        fi
        printf '%s\\n' "\$cleaned" >> cleaned_proteins.list
    done < protein_fastas.list
    cleaned_proteins=\$(paste -sd, cleaned_proteins.list)

    braker_cmd=(
        "\${BRAKER_DIR}/scripts/braker.pl"
        "--genome=${genome}"
        "--bam=\${bam}"
        "--prot_seq=\${cleaned_proteins}"
        "--threads=${task.cpus}"
        "--workingdir=\${PWD}"
        "--softmasking"
        "--gff3"
        "--PROTHINT_PATH=\${PROTHINT_BIN}/"
        "--GENEMARK_PATH=\${GENEMARK_DIR}"
        "--AUGUSTUS_CONFIG_PATH=\${AUGUSTUS_CONFIG_DIR}"
        "--TSEBRA_PATH=\${TSEBRA_BIN}"
    )

    printf '%q ' "\${braker_cmd[@]}" > braker3_command.txt
    printf '\\n' >> braker3_command.txt
    log "Executing BRAKER3 command recorded in braker3_command.txt"

    "\${braker_cmd[@]}"
    # to test to decrease monoexon genes : --augustus_args="--singlestrand=true --alternatives-from-evidence=0"
    # --singlestrand=true: If enabled, only keeps predictions consistent with a single strand.
    # --alternatives-from-evidence=0: Disables alternative predictions based on hints, which can help filter out monoexons.
    copy_required_output "Augustus/augustus.hints.gff3" "augustus.hints.gff3"
    copy_required_output "GeneMark-ETP/genemark.gtf" "genemark.gtf"
    copy_required_output "GeneMark-ETP/genemark_supported.gtf" "genemark_supported.gtf"
    copy_required_output "braker.gff3" "braker.gff3"

    {
        printf '"%s":\\n' "${task.process}"
        printf '  braker: "3.0.8"\\n'
        printf '  prothint: "2.6.0"\\n'
        printf '  genemark_etp_path: "%s"\\n' "\${GENEMARK_DIR}"
        printf '  augustus_config_path: "%s"\\n' "\${AUGUSTUS_CONFIG_DIR}"
        printf '  tsebra_path: "%s"\\n' "\${TSEBRA_BIN}"
        printf '  bam_short_count: "%s"\\n' "\${bam_short_count}"
        printf '  bam_long_count: "%s"\\n' "\${bam_long_count}"
        printf '  bam_all_count: "%s"\\n' "\${bam_all_count}"
        printf '  protein_count: "%s"\\n' "\${protein_count}"
    } > versions.yml
    log "BRAKER3 completed successfully"
    """

  stub:
    """
    printf "##gff-version 3\\n" > augustus.hints.gff3
    : > genemark.gtf
    cp genemark.gtf genemark_supported.gtf
    printf "##gff-version 3\\n" > braker.gff3
    printf "Stub BRAKER3 run with long reads\\n" > braker3_run.log
    printf "stub\\n" > braker3_command.txt
    printf "type\\tpath\\tsize_bytes\\n" > braker3_inputs.tsv
    printf '"%s":\n  braker: "stub"\n' "${task.process}" > versions.yml
    """
}
