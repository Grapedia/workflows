process EDTA {

  tag "Executing EDTA TE annotation on the following genome: $genome"
  container 'avelt/edta:latest'
  containerOptions "--volume ${genome_path}:/genome_path --volume ${projectDir}/scripts/:/scripts"
  publishDir "$projectDir/FINAL_OUTPUT"
  cpus 5

  input:
    val(genome_path)
    val(genome)

  output:
    path("*TElib.fa"), emit : TElib_fasta
    path("*TEanno.gff3"), emit : TE_annotations_gff3
    path("*MAKER.masked"), emit : masked_genome

  script:
    """
    /scripts/edta.sh -g /genome_path/$genome -n ${task.cpus}
    """
}
