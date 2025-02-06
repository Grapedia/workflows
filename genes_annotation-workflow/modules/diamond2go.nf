process diamond2go {

  tag "Executing diamond2go annotation on $proteins_file"
  container 'avelt/diamond2go:latest'
  publishDir "$projectDir/FINAL_OUTPUT/diamond2go/"
  cpus 5

  input:
    path(proteins_file)

  output:
    path("*-diamond*")

  script:
    """
    DATE=\$(date "+%Y-%m-%d %H:%M:%S")
    echo "[\$DATE] Running diamond2go on $proteins_file"
    CMD="perl /Diamond2GO/Diamond2go.pl -d /Diamond2GO/resources/nr_clean_d2go.dmnd -q $proteins_file -t protein"
    echo "[\$DATE] Executing: \$CMD"
    perl /Diamond2GO/Diamond2go.pl -d /Diamond2GO/resources/nr_clean_d2go.dmnd -q $proteins_file -t protein
    """
}
