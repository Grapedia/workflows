// 1. Generate reference genome indices using HISAT2
process hisat2_genome_indices {
  label 'process_index'

  tag "HISAT2 indexes generation on $genome"
  container params.container_hisat2
  publishDir "${params.output_dir}/intermediate_files/evidence_data/hisat2_databases/", mode: "copy", enabled: params.publish_intermediates
  input:
    path(genome_fasta)
    val(genome)

  output:
    path("hisat2_index"), emit: index
    path "versions.yml", emit: versions



  script:
    """
    set -euo pipefail
    DATE=\$(date "+%Y-%m-%d %H:%M:%S")
    echo "[\$DATE] Running HISAT2 index creation on $genome"
    mkdir -p hisat2_index
    CMD="/hisat2-2.2.1/hisat2-build -p ${task.cpus} ${genome_fasta} hisat2_index/$genome"
    echo "[\$DATE] Executing: \$CMD"
    /hisat2-2.2.1/hisat2-build -p ${task.cpus} ${genome_fasta} hisat2_index/$genome
    /hisat2-2.2.1/hisat2-build --version 2>&1 | head -n 1 | sed 's/^/  hisat2_build: "/; s/\$/"/' | {
      printf '"%s":\\n' "${task.process}"
      cat
    } > versions.yml
    """

  stub:
    """
    set -euo pipefail

    mkdir -p hisat2_index
    for i in 1 2 3 4 5 6 7 8; do
      printf "stub HISAT2 index\\n" > hisat2_index/${genome}.\${i}.ht2
    done
    printf '"%s":\n  hisat2_build: "stub"\n' "${task.process}" > versions.yml
    """
}
