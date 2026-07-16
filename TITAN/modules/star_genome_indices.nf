// 1. Generate reference genome indices using STAR
process star_genome_indices {
  label 'process_index'

  tag "STAR genomeGenerate on $genome"
  container params.container_star
  publishDir "${params.output_dir}/intermediate_files/evidence_data/star_databases/", mode: "copy", enabled: params.publish_intermediates
  input:
    path(genome_fasta)
    val(genome)

  output:
    path("${genome}_index"), emit: index
    path "versions.yml", emit: versions



  script:
    """
    set -euo pipefail
    DATE=\$(date "+%Y-%m-%d %H:%M:%S")
    echo "[\$DATE] Running STAR index creation on $genome"
    CMD="STAR --runThreadN ${task.cpus} --runMode genomeGenerate --genomeDir ${genome}_index --genomeFastaFiles ${genome_fasta}"
    echo "[\$DATE] Executing: \$CMD"
    STAR --runThreadN ${task.cpus} --runMode genomeGenerate --genomeDir ${genome}_index --genomeFastaFiles ${genome_fasta}
    chmod -R 755 ${genome}_index
    printf '"%s":\n  container: "not_recorded"\n' "${task.process}" > versions.yml
    """

  stub:
    """
    mkdir -p ${genome}_index
    printf "stub STAR index\\n" > ${genome}_index/Genome
    printf '"%s":\n  container: "not_recorded"\n' "${task.process}" > versions.yml
    """
}
