// 2. Generate minimap2 genome indices
process minimap2_genome_indices {

  tag "Minimap2 genome indices on ${genome}"
  container params.container_minimap2_samtools
  containerOptions "--volume $genome_path:/genome_path"
  publishDir "${params.output_dir}/intermediate_files/evidence_data/minimap2_databases/"
  cpus 4

  input:
    val(genome_path)
    val(genome)

  output:
    path("${genome}.mmi"), emit : index

  script:
    """
    DATE=\$(date "+%Y-%m-%d %H:%M:%S")
    echo "[\$DATE] Running minimap2 index creation on $genome"
    CMD="minimap2 -d ${genome}.mmi /genome_path/$genome"
    echo "[\$DATE] Executing: \$CMD"
    minimap2 -d ${genome}.mmi /genome_path/$genome
    """

  stub:
    """
    printf "stub minimap2 index\\n" > ${genome}.mmi
    """
}
