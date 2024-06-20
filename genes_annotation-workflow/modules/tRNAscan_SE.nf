process tRNAscan_SE {

  tag "Executing tRNAscan-SE annotation on the following genome: $genome"
  container 'quay.io/biocontainers/trnascan-se:2.0--pl526h470a237_1'
  containerOptions "--volume ${genome_path}:/genome_path"
  publishDir "$projectDir/FINAL_OUTPUT"
  cpus 4

  input:
    val(genome_path)
    val(genome)

  output:
    file("annotations.tRNAscan-SE.txt")
    file("summary.txt")

  script:
    """
    tRNAscan-SE -o annotations.tRNAscan-SE.txt -m summary.txt -d -E -p tRNAscan-SE /genome_path/$genome
    """
}
