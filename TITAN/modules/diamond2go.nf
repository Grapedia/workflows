process diamond2go {

  tag "Executing diamond2go annotation on $proteins_file"
  container 'avelt/diamond2go:latest'
  publishDir "${params.output_dir}/Diamond2GO_outputs", mode: 'copy'
  cpus 5

  input:
    path(proteins_file_all)
    path(proteins_file_main)

  output:
    path("*-diamond*")

  script:
    """
    DATE=\$(date "+%Y-%m-%d %H:%M:%S")
    
    echo "[\$DATE] Running diamond2go on $proteins_file_all"
    CMD="perl /Diamond2GO/Diamond2go.pl -d /Diamond2GO/resources/nr_clean_d2go.dmnd -q $proteins_file_all -t protein"
    echo "[\$DATE] Executing: \$CMD"
    perl /Diamond2GO/Diamond2go.pl -d /Diamond2GO/resources/nr_clean_d2go.dmnd -q $proteins_file_all -t protein

    echo "[\$DATE] Running diamond2go on $proteins_file_main"
    CMD="perl /Diamond2GO/Diamond2go.pl -d /Diamond2GO/resources/nr_clean_d2go.dmnd -q $proteins_file_main -t protein"
    echo "[\$DATE] Executing: \$CMD"
    perl /Diamond2GO/Diamond2go.pl -d /Diamond2GO/resources/nr_clean_d2go.dmnd -q $proteins_file_main -t protein
    """
}
