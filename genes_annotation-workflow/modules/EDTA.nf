process EDTA {

  tag "Executing EDTA TE annotation on the following genome: $genome"
  container 'avelt/edta:latest'
  containerOptions "--volume ${genome_path}:/genome_path --volume ${projectDir}/scripts/:/scripts"
  publishDir "$projectDir/FINAL_OUTPUT"
  cpus 10

  input:
    val(genome_path)
    val(genome)

  output:
    file("*.fa"), emit : TElib_fasta
    file("*.gff3"), emit : TE_annotations_gff3

  script:
    """
    /scripts/edta.sh -g /genome_path/$genome -n ${task.cpus}
    """
}
