process EDTA {
  label 'process_prediction'

  tag "Executing EDTA TE annotation on the following genome: $genome"
  container params.container_edta
  publishDir "${params.output_dir}/tmp", mode: 'copy', enabled: params.publish_intermediates
  publishDir "${params.output_dir}", mode: 'copy', saveAs: { filename ->
    filename.endsWith('.MAKER.masked') ? 'assembly_masked.EDTA.fasta' : null
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
    cp ./*TElib.fa edta.TElib.fa
    cp ./*TEanno.gff3 edta.TEanno.gff3
    cp ./*MAKER.masked edta.MAKER.masked
    printf '"%s":\n  container: "not_recorded"\n' "${task.process}" > versions.yml
    """

  stub:
    """
    cp ${genome_fasta} edta.MAKER.masked
    printf ">stub_te\\nNNNN\\n" > edta.TElib.fa
    printf "##gff-version 3\\n" > edta.TEanno.gff3
    printf '"%s":\n  container: "not_recorded"\n' "${task.process}" > versions.yml
    """
}
