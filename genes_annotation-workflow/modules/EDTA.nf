process EDTA {

  tag "Executing EDTA TE annotation on the following genome: $genome"
  container 'quay.io/biocontainers/edta:2.2.0--hdfd78af_1'
  containerOptions "--volume ${genome_path}:/genome_path"
  publishDir "$projectDir/FINAL_OUTPUT"
  cpus 4

  input:
    val(genome_path)
    val(genome)

  output:
    file("annotations.EDTA.txt")

  script:
    """
    perl /usr/local/bin/EDTA.pl --genome /genome_path/$genome --species others --step all --sensitive 1 --anno 1 --threads ${task.cpus}
    """
}
