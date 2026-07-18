// Gene-set completeness of the final AEGIS annotation, protein mode.
process busco {
  label 'process_aegis'

  tag "Executing BUSCO on final AEGIS main protein set"
  container params.container_busco
  publishDir "${params.output_dir}/quality_report/busco", mode: 'copy', saveAs: { filename ->
    if (filename in [
      'busco_short_summary.txt',
      'busco_full_table.tsv',
      'busco_missing_busco_list.tsv',
      'versions.yml'
    ]) {
      return filename
    }
    return null
  }

  input:
    path(proteins_file_main)

  output:
    path "busco_short_summary.txt", emit: short_summary
    path "busco_full_table.tsv", emit: full_table
    path "busco_missing_busco_list.tsv", emit: missing_list
    path "versions.yml", emit: versions



  script:
    """
    set -euo pipefail

    if [[ "${params.run_busco}" != "true" ]]; then
        touch busco_short_summary.txt
        touch busco_full_table.tsv
        touch busco_missing_busco_list.tsv
        printf '"%s":\n  busco: "skipped"\n  container: "%s"\n' \
            "${task.process}" "${task.container}" > versions.yml
        exit 0
    fi

    if [[ ! -s "${proteins_file_main}" ]]; then
        echo "ERROR: missing or empty main proteins FASTA: ${proteins_file_main}" >&2
        exit 1
    fi

    if [[ -z "${params.busco_data_dir}" || "${params.busco_data_dir}" == "false" ]]; then
        echo "ERROR: busco_data_dir must be provided when run_busco=true" >&2
        exit 1
    fi

    busco -i "${proteins_file_main}" -m protein -l "${params.busco_lineage}" \\
      -o busco_out -c ${task.cpus} --offline --download_path "${params.busco_data_dir}" -f

    cp busco_out/short_summary.*.busco_out.txt busco_short_summary.txt
    cp "busco_out/run_${params.busco_lineage}/full_table.tsv" busco_full_table.tsv
    cp "busco_out/run_${params.busco_lineage}/missing_busco_list.tsv" busco_missing_busco_list.tsv

    test -s busco_short_summary.txt

    printf '"%s":\n  busco: "container-pinned"\n  lineage: "%s"\n  container: "%s"\n' \\
        "${task.process}" "${params.busco_lineage}" "${task.container}" > versions.yml
    """

  stub:
    """
    set -euo pipefail
    printf "# BUSCO version is: stub\\nC:0%%\\n" > busco_short_summary.txt
    printf "Busco_id\\tStatus\\n" > busco_full_table.tsv
    printf "Busco_id\\n" > busco_missing_busco_list.tsv
    printf '"%s":\n  busco: "stub"\n' "${task.process}" > versions.yml
    """
}
