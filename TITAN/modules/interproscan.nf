process interproscan {
  label 'process_aegis'

  tag "Executing InterProScan on $proteins_file_all and $proteins_file_main"
  container params.container_interproscan
  publishDir "${params.output_dir}/InterProScan_outputs", mode: 'copy', saveAs: { filename ->
    if (filename in [
      'final_annotation_proteins_all.tsv',
      'final_annotation_proteins_main.tsv',
      'final_annotation_proteins_all.gff3',
      'final_annotation_proteins_main.gff3',
      'final_annotation_proteins_all.json',
      'final_annotation_proteins_main.json',
      'versions.yml'
    ]) {
      return filename
    }
    return null
  }

  input:
    path(proteins_file_all)
    path(proteins_file_main)

  output:
    path "final_annotation_proteins_all.tsv", emit: proteins_all_tsv
    path "final_annotation_proteins_main.tsv", emit: proteins_main_tsv
    path "final_annotation_proteins_all.gff3", emit: proteins_all_gff3
    path "final_annotation_proteins_main.gff3", emit: proteins_main_gff3
    path "final_annotation_proteins_all.json", emit: proteins_all_json
    path "final_annotation_proteins_main.json", emit: proteins_main_json
    path "versions.yml", emit: versions

  script:
    """
    set -euo pipefail

    IPRSCAN="/opt/interproscan/interproscan.sh"

    require_file() {
        local label="\$1"
        local file="\$2"
        if [[ ! -s "\$file" ]]; then
            echo "ERROR: missing or empty \${label}: \${file}" >&2
            exit 1
        fi
    }

    if [[ "${params.run_interproscan}" != "true" ]]; then
        touch final_annotation_proteins_all.tsv
        touch final_annotation_proteins_main.tsv
        touch final_annotation_proteins_all.gff3
        touch final_annotation_proteins_main.gff3
        touch final_annotation_proteins_all.json
        touch final_annotation_proteins_main.json
        printf '"%s":\n  interproscan: "skipped"\n  container: "%s"\n' \
            "${task.process}" "${task.container}" > versions.yml
        exit 0
    fi

    require_file "all proteins FASTA" "${proteins_file_all}"
    require_file "main proteins FASTA" "${proteins_file_main}"

    if [[ -z "${params.interproscan_data_dir}" || "${params.interproscan_data_dir}" == "false" ]]; then
        echo "ERROR: interproscan_data_dir must be provided when run_interproscan=true" >&2
        exit 1
    fi

    mkdir -p interproscan_tmp

    "\$IPRSCAN" -i "${proteins_file_all}" -b final_annotation_proteins_all -T interproscan_tmp -cpu ${task.cpus} -dp -f TSV,GFF3,JSON -goterms -pathways
    "\$IPRSCAN" -i "${proteins_file_main}" -b final_annotation_proteins_main -T interproscan_tmp -cpu ${task.cpus} -dp -f TSV,GFF3,JSON -goterms -pathways

    # TSV can legitimately be empty when a protein has zero domain hits across
    # all analyses; GFF3 always carries header lines, so check that instead.
    test -s final_annotation_proteins_all.gff3
    test -s final_annotation_proteins_main.gff3

    interproscan_version=\$("\$IPRSCAN" -version 2>/dev/null | head -n 1 || true)
    printf '"%s":\n  interproscan: "%s"\n  container: "%s"\n' \
        "${task.process}" "\${interproscan_version:-unknown}" "${task.container}" > versions.yml
    """

  stub:
    """
    set -euo pipefail
    printf "query\\tid\\n" > final_annotation_proteins_all.tsv
    printf "query\\tid\\n" > final_annotation_proteins_main.tsv
    printf "##gff-version 3\\n" > final_annotation_proteins_all.gff3
    printf "##gff-version 3\\n" > final_annotation_proteins_main.gff3
    printf "[]\\n" > final_annotation_proteins_all.json
    printf "[]\\n" > final_annotation_proteins_main.json
    printf '"%s":\n  interproscan: "stub"\n  container: "%s"\n' "${task.process}" "${task.container}" > versions.yml
    """
}
