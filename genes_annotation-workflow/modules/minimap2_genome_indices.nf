// 2. Generate minimap2 genome indices
process minimap2_genome_indices {

  tag "Minimap2 genome indices on ${genome}"
  container 'avelt/minimap2_samtools:latest'
  containerOptions "--volume $genome_path:/genome_path"
  publishDir "$params.outdir/evidence_data/minimap2_databases/"
  cpus 4

  input:
    val(genome_path)
    val(genome)

  output:
    file("${genome}.mmi")

  when:
  params.use_long_reads

  script:
    """
    DATE=\$(date "+%Y-%m-%d %H:%M:%S")
    echo "[\$DATE] Running minimap2 index creation on $genome" >> ${params.logfile} 2>&1
    CMD="minimap2 -d ${genome}.mmi /genome_path/$genome"
    echo "[\$DATE] Executing: \$CMD" >> ${params.logfile} 2>&1
    minimap2 -d ${genome}.mmi /genome_path/$genome
    """
}
