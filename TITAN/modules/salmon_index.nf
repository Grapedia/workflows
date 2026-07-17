process salmon_index {

  label 'process_index'
  tag "Executing salmon indexing on $cds_fasta"
  container params.container_salmon
  input:
    path(cds_fasta)

  output:
    path("salmon_index"), emit: index
    path "versions.yml", emit: versions



  script:
    """
    set -euo pipefail
    DATE=\$(date "+%Y-%m-%d %H:%M:%S")
    echo "[\$DATE] Running salmon indexing on ${cds_fasta}"
    CMD="salmon index -t ${cds_fasta} -i salmon_index"
    echo "[\$DATE] Executing: \$CMD"
    salmon index -t ${cds_fasta} -i salmon_index
    salmon --version | sed 's/^/  salmon: "/; s/\$/"/' | {
      printf '"%s":\\n' "${task.process}"
      cat
    } > versions.yml
    """

  stub:
    """
    set -euo pipefail

    mkdir -p salmon_index
    printf "stub salmon index\\n" > salmon_index/versionInfo.json
    printf '"%s":\n  salmon: "stub"\n' "${task.process}" > versions.yml
    """
}
