// Rfam covariance-model ncRNA annotation on the target genome with Infernal.
process infernal_rfam {
  label 'process_rfam'

  tag "Infernal/Rfam ncRNA annotation"
  container params.container_infernal
  publishDir "${params.output_dir}/additional_annotations/ncrna/rfam", mode: 'copy'

  input:
    path(genome)
    path(rfam_tblout_to_gff3)

  output:
    path "rfam_hits.tbl", emit: tblout
    path "rfam_search.out", emit: search_log
    path "rfam_ncrna.gff3", emit: gff3
    path "versions.yml", emit: versions

  script:
    """
    set -euo pipefail

    if [[ "${params.run_rfam}" != "true" ]]; then
        touch rfam_hits.tbl rfam_search.out
        printf "##gff-version 3\\n" > rfam_ncrna.gff3
        printf '"%s":\n  infernal: "skipped"\n  rfam_data_dir: "%s"\n  container: "%s"\n' \\
            "${task.process}" "${params.rfam_data_dir}" "${task.container}" > versions.yml
        exit 0
    fi

    if [[ ! -s "${genome}" ]]; then
        echo "ERROR: missing or empty target genome FASTA: ${genome}" >&2
        exit 1
    fi

    if [[ -z "${params.rfam_data_dir}" || "${params.rfam_data_dir}" == "false" ]]; then
        echo "ERROR: rfam_data_dir must be provided when run_rfam=true" >&2
        exit 1
    fi

    if [[ ! -s "${params.rfam_data_dir}/Rfam.cm" || ! -s "${params.rfam_data_dir}/Rfam.clanin" ]]; then
        echo "ERROR: rfam_data_dir must contain Rfam.cm and Rfam.clanin" >&2
        exit 1
    fi

    cmsearch --cpu ${task.cpus} \\
      --tblout rfam_hits.tbl \\
      --fmt 2 \\
      --cut_ga \\
      --rfam \\
      --nohmmonly \\
      --clanin "${params.rfam_data_dir}/Rfam.clanin" \\
      "${params.rfam_data_dir}/Rfam.cm" \\
      "${genome}" > rfam_search.out

    python3 "${rfam_tblout_to_gff3}" rfam_hits.tbl > rfam_ncrna.gff3
    test -s rfam_ncrna.gff3

    cmsearch_version=\$(cmsearch -h 2>&1 | head -n 2 | tr '\\n' ' ' | sed 's/"/\\\\\\"/g')
    rfam_release="unknown"
    if [[ -s "${params.rfam_data_dir}/Rfam.version" ]]; then
        rfam_release=\$(tr '\\n' ' ' < "${params.rfam_data_dir}/Rfam.version" | sed 's/"/\\\\\\"/g')
    fi
    printf '"%s":\n  infernal: "%s"\n  rfam_data_dir: "%s"\n  rfam_release: "%s"\n  container: "%s"\n' \\
        "${task.process}" "\${cmsearch_version}" "${params.rfam_data_dir}" "\${rfam_release}" "${task.container}" > versions.yml
    """

  stub:
    """
    set -euo pipefail
    cat > rfam_hits.tbl <<'EOF'
#target name        accession  query name  accession  mdl  mdl from  mdl to  seq from  seq to  strand  trunc  pass  gc  bias  score  E-value  inc  description of target
chrStub             -          5S_rRNA     RF00001    cm   1         119     20        138     +       no     1     0.49 0.0   87.4   1e-20    !    stub 5S ribosomal RNA
EOF
    printf "# Infernal/Rfam stub\\n" > rfam_search.out
    printf "##gff-version 3\\nchrStub\\tInfernal/Rfam\\trRNA\\t20\\t138\\t87.4\\t+\\t.\\tID=rfam.chrStub.RF00001.1;Name=5S_rRNA;Rfam_ID=RF00001;E_value=1e-20;rfam_type=rRNA\\n" > rfam_ncrna.gff3
    printf '"%s":\n  infernal: "stub"\n' "${task.process}" > versions.yml
    """
}
