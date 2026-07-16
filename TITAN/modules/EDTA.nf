process EDTA {
  label 'process_prediction'

  tag "Executing EDTA TE annotation on the following genome: $genome"
  container params.container_edta
  publishDir "${params.output_dir}/tmp", mode: 'copy'
  publishDir "${params.output_dir}", mode: 'copy', saveAs: { filename ->
    filename.endsWith('.MAKER.masked') ? 'assembly_masked.EDTA.fasta' : null
  }

  input:
    path(genome_fasta)
    val(genome)
    path(edta_script)

  output:
    path("*TElib.fa"), emit : TElib_fasta
    path("*TEanno.gff3"), emit : TE_annotations_gff3
    path("*MAKER.masked"), emit : masked_genome

  script:
    """
    DATE=\$(date "+%Y-%m-%d %H:%M:%S")
    echo "[\$DATE] Running EDTA on $genome"
    CMD="${edta_script} -g ${genome_fasta} -n ${task.cpus}"
    echo "[\$DATE] Executing: \$CMD"
    ${edta_script} -g ${genome_fasta} -n ${task.cpus}
    """

  stub:
    """
    printf ">stub_masked\\nNNNN\\n" > ${genome}.MAKER.masked
    printf ">stub_te\\nNNNN\\n" > ${genome}.TElib.fa
    printf "##gff-version 3\\n" > ${genome}.TEanno.gff3
    """
}
