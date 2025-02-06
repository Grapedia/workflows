// 1. Generate reference genome indices using HISAT2
process hisat2_genome_indices {

  tag "HISAT2 indexes generation on $genome"
  container 'avelt/hisat2:latest'
  containerOptions "--volume $genome_path:/genome_path"
  publishDir "$params.outdir/evidence_data/hisat2_databases/"
  cpus 4

  input:
    val(genome_path)
    val(genome)

  output:
    path("${genome}.*.ht2")

  script:
    """
    DATE=\$(date "+%Y-%m-%d %H:%M:%S")
    echo "[\$DATE] Running HISAT2 index creation on $genome"
    CMD="/hisat2-2.2.1/hisat2-build -p ${task.cpus} /genome_path/$genome $genome"
    echo "[\$DATE] Executing: \$CMD"
    /hisat2-2.2.1/hisat2-build -p ${task.cpus} /genome_path/$genome $genome
    chmod -R 755 ${genome}*
    """
}
