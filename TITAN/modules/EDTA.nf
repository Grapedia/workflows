process EDTA {
  label 'process_prediction'

  tag "Executing EDTA TE annotation on the following genome: $genome"
  container params.container_edta
  publishDir "${params.output_dir}", mode: 'copy', saveAs: { filename ->
    if (filename == 'edta.MAKER.masked') {
      return 'assembly_masked.EDTA.fasta'
    }
    if (params.publish_intermediates && filename in ['edta.TElib.fa', 'edta.TEanno.gff3', 'versions.yml']) {
      return "intermediate_files/evidence_data/EDTA/${filename}"
    }
    return null
  }

  input:
    path(genome_fasta)
    val(genome)
    path(edta_script)

  output:
    path("edta.TElib.fa"), emit: TElib_fasta
    path("edta.TEanno.gff3"), emit: TE_annotations_gff3
    path("edta.MAKER.masked"), emit: masked_genome
    path "versions.yml", emit: versions



  script:
    """
    set -euo pipefail
    DATE=\$(date "+%Y-%m-%d %H:%M:%S")
    echo "[\$DATE] Running EDTA on $genome"
    CMD="${edta_script} -g ${genome_fasta} -n ${task.cpus}"
    echo "[\$DATE] Executing: \$CMD"
    ${edta_script} -g ${genome_fasta} -n ${task.cpus}

    copy_unique() {
      local pattern="\$1"
      local output="\$2"
      mapfile -t matches < <(find . -maxdepth 1 -type f -name "\$pattern" | sort)
      if [[ "\${#matches[@]}" -ne 1 ]]; then
        printf 'ERROR: expected exactly one EDTA output matching %s, found %s\\n' "\$pattern" "\${#matches[@]}" >&2
        printf '  %s\\n' "\${matches[@]}" >&2
        exit 1
      fi
      cp "\${matches[0]}" "\$output"
    }

    copy_unique "*TElib.fa" edta.TElib.fa
    copy_unique "*TEanno.gff3" edta.TEanno.gff3
    copy_unique "*MAKER.masked" edta.MAKER.masked
    printf '"%s":\n  container: "not_recorded"\n' "${task.process}" > versions.yml
    """

  stub:
    """
    set -euo pipefail
    cp ${genome_fasta} edta.MAKER.masked
    printf ">stub_te\\nNNNN\\n" > edta.TElib.fa
    printf "##gff-version 3\\n" > edta.TEanno.gff3
    printf '"%s":\n  container: "not_recorded"\n' "${task.process}" > versions.yml
    """
}
