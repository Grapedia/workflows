process salmon_index {

  label 'process_index'
  tag "Executing salmon indexing on $cds_fasta"
  container params.container_salmon
  publishDir "${params.output_dir}/intermediate_files/salmon_index/"
  input:
    path(cds_fasta)

  output:
    path("salmon_index"), emit : index
    path "versions.yml", emit: versions



  script:
    """
    set -euo pipefail
    DATE=\$(date "+%Y-%m-%d %H:%M:%S")
    echo "[\$DATE] Running salmon indexing on ${cds_fasta}"
    CMD="salmon index -t ${cds_fasta} -i salmon_index"
    echo "[\$DATE] Executing: \$CMD"
    salmon index -t ${cds_fasta} -i salmon_index
    printf '"%s":\n  container: "not_recorded"\n' "${task.process}" > versions.yml
    """

  stub:
    """
    mkdir -p salmon_index
    printf "stub salmon index\\n" > salmon_index/versionInfo.json
    printf '"%s":\n  container: "not_recorded"\n' "${task.process}" > versions.yml
    """
}
