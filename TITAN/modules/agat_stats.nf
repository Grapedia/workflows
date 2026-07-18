// Structural statistics (gene/mRNA/exon counts and lengths) of the final
// AEGIS annotation, using the same AGAT image already pinned for CDS
// extraction elsewhere in the pipeline.
process agat_stats {
  label 'process_low'

  tag "AGAT structural statistics on the final AEGIS annotation"
  container params.container_agat
  publishDir "${params.output_dir}/quality_report/agat_stats", mode: 'copy'

  input:
    path(annotation_gff3)

  output:
    path "agat_stats.txt", emit: stats_txt
    path "versions.yml", emit: versions



  script:
    """
    set -euo pipefail
    agat_sp_statistics.pl --gff ${annotation_gff3} --output agat_stats.txt

    printf '"%s":\n  agat_sp_statistics: "container-pinned"\n  container: "%s"\n' \\
        "${task.process}" "${task.container}" > versions.yml
    """

  stub:
    """
    set -euo pipefail
    printf "Number of gene\\t0\\n" > agat_stats.txt
    printf '"%s":\n  agat_sp_statistics: "stub"\n' "${task.process}" > versions.yml
    """
}
