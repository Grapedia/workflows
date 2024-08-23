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

  script:
    """
    minimap2 -d ${genome}.mmi /genome_path/$genome
    """
}
