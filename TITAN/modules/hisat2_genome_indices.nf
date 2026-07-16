// 1. Generate reference genome indices using HISAT2
process hisat2_genome_indices {
  label 'process_index'

  tag "HISAT2 indexes generation on $genome"
  container params.container_hisat2
  publishDir "${params.output_dir}/intermediate_files/evidence_data/hisat2_databases/"
  input:
    path(genome_fasta)
    val(genome)

  output:
    path("${genome}.*.ht2"), emit : index
    path "versions.yml", emit: versions



  script:
    """
    set -euo pipefail
    DATE=\$(date "+%Y-%m-%d %H:%M:%S")
    echo "[\$DATE] Running HISAT2 index creation on $genome"
    CMD="/hisat2-2.2.1/hisat2-build -p ${task.cpus} ${genome_fasta} $genome"
    echo "[\$DATE] Executing: \$CMD"
    /hisat2-2.2.1/hisat2-build -p ${task.cpus} ${genome_fasta} $genome
    chmod -R 755 ${genome}*
    printf '"%s":\n  container: "not_recorded"\n' "${task.process}" > versions.yml
    """

  stub:
    """
    for i in 1 2 3 4 5 6 7 8; do
      printf "stub HISAT2 index\\n" > ${genome}.\${i}.ht2
    done
    printf '"%s":\n  container: "not_recorded"\n' "${task.process}" > versions.yml
    """
}
