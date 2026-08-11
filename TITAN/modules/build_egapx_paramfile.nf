process build_egapx_paramfile {
  label 'process_low'
  tag "EGAPx paramfile from fastp-trimmed short reads"
  publishDir "${params.output_dir}/intermediate_files/evidence_data/egapx", mode: "copy", enabled: params.publish_intermediates

  input:
    path egapx_paramfile
    val trimmed_rows
    path trimmed_fastqs

  output:
    path "egapx.trimmed_reads.yaml", emit: paramfile
    path "egapx.trimmed_reads.tsv", emit: manifest
    path "versions.yml", emit: versions

  script:
    def manifestRows = trimmed_rows.collect { row ->
      def sample_ID = row[0]
      def library_layout = row[1]
      if (library_layout == 'paired') {
        return "${sample_ID}\t${library_layout}\t${sample_ID}_1.trimmed.fastq.gz\t${sample_ID}_2.trimmed.fastq.gz"
      }
      return "${sample_ID}\t${library_layout}\t${sample_ID}_1.trimmed.fastq.gz\t"
    }.join('\n')

    """
    set -euo pipefail

    cat > egapx.trimmed_reads.tsv <<'EOF'
${manifestRows}
EOF

    python3 ${projectDir}/scripts/build_egapx_paramfile.py \\
      --input "${egapx_paramfile}" \\
      --trimmed-reads egapx.trimmed_reads.tsv \\
      --output egapx.trimmed_reads.yaml

    printf '"%s":\n  egapx_paramfile_source: "%s"\n  short_reads_source: "fastp_trimmed"\n' \\
      "${task.process}" "${egapx_paramfile}" > versions.yml
    """

  stub:
    def manifestRows = trimmed_rows.collect { row ->
      def sample_ID = row[0]
      def library_layout = row[1]
      if (library_layout == 'paired') {
        return "${sample_ID}\t${library_layout}\t${sample_ID}_1.trimmed.fastq.gz\t${sample_ID}_2.trimmed.fastq.gz"
      }
      return "${sample_ID}\t${library_layout}\t${sample_ID}_1.trimmed.fastq.gz\t"
    }.join('\n')

    """
    set -euo pipefail

    cat > egapx.trimmed_reads.tsv <<'EOF'
${manifestRows}
EOF

    python3 ${projectDir}/scripts/build_egapx_paramfile.py \\
      --input "${egapx_paramfile}" \\
      --trimmed-reads egapx.trimmed_reads.tsv \\
      --output egapx.trimmed_reads.yaml

    printf '"%s":\n  egapx_paramfile_source: "%s"\n  short_reads_source: "fastp_trimmed"\n' \\
      "${task.process}" "${egapx_paramfile}" > versions.yml
    """
}
