process EDTA {

  tag "Executing EDTA TE annotation on the following genome: $genome"
  container params.container_edta
  containerOptions "--memory=50g --volume ${genome_path}:/genome_path --volume ${projectDir}/scripts/:/scripts"
  cpus params.edta_cpus
  publishDir "${params.output_dir}/tmp", mode: 'copy'
  publishDir "${params.output_dir}", mode: 'copy', saveAs: { filename ->
    filename.endsWith('.MAKER.masked') ? 'assembly_masked.EDTA.fasta' : null
  }

  input:
    val(genome_path)
    val(genome)

  output:
    path("*TElib.fa"), emit : TElib_fasta
    path("*TEanno.gff3"), emit : TE_annotations_gff3
    path("*MAKER.masked"), emit : masked_genome

  script:
    """
    DATE=\$(date "+%Y-%m-%d %H:%M:%S")
    echo "[\$DATE] Running EDTA on $genome"
    CMD="/scripts/edta.sh -g /genome_path/$genome -n ${task.cpus}"
    echo "[\$DATE] Executing: \$CMD"
    /scripts/edta.sh -g /genome_path/$genome -n ${task.cpus}
    """

  stub:
    """
    printf ">stub_masked\\nNNNN\\n" > ${genome}.MAKER.masked
    printf ">stub_te\\nNNNN\\n" > ${genome}.TElib.fa
    printf "##gff-version 3\\n" > ${genome}.TEanno.gff3
    """
}
