process diamond2go {

  tag "Executing diamond2go annotation on $proteins_file"
  container 'avelt/diamond2go:latest'
  publishDir "$projectDir/FINAL_OUTPUT/diamond2go/"
  cpus 5

  input:
    val(proteins_file)

  output:
    path("*-diamond*")

  script:
    """
    perl /Diamond2GO/Diamond2go.pl -d /Diamond2GO/resources/nr_clean_d2go.dmnd -q $proteins_file -t protein
    """
}
