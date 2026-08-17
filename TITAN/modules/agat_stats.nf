// Structural statistics for the final AEGIS annotation.
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

// Same AGAT stats as above, on the "high confidence" GFF3 (unsupported
// single-exon genes removed by flag_low_confidence_monoexonic_genes) so its
// gene-set structure can be compared side by side against the full
// annotation's. A distinct process (not a parameterized reuse of
// `agat_stats`) so this never changes the task hash / cache of the
// already-validated primary agat_stats run.
process agat_stats_high_confidence_monoexonic {
  label 'process_low'

  tag "AGAT structural statistics on the high-confidence (monoexonic-filtered) annotation"
  container params.container_agat
  publishDir "${params.output_dir}/aegis_outputs_high_confidence_monoexonic/agat_stats", mode: 'copy'

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
