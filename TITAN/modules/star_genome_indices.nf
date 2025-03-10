// 1. Generate reference genome indices using STAR
process star_genome_indices {

  tag "STAR genomeGenerate on $genome"
  container 'quay.io/biocontainers/star:2.7.11b--h43eeafb_2'
  containerOptions "--volume $genome_path:/genome_path"
  publishDir "${params.output_dir}/intermediate_files/evidence_data/star_databases/"
  cpus 4

  input:
    val(genome_path)
    val(genome)

  output:
    path("${genome}_index"), emit: index

  script:
    """
    DATE=\$(date "+%Y-%m-%d %H:%M:%S")
    echo "[\$DATE] Running STAR index creation on $genome"
    CMD="STAR --runThreadN ${task.cpus} --runMode genomeGenerate --genomeDir ${genome}_index --genomeFastaFiles /genome_path/$genome"
    echo "[\$DATE] Executing: \$CMD"
    STAR --runThreadN ${task.cpus} --runMode genomeGenerate --genomeDir ${genome}_index --genomeFastaFiles /genome_path/$genome
    chmod -R 755 ${genome}_index
    """
}
