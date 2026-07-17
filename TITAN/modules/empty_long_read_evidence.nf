process empty_long_read_evidence {
  label 'process_low'
  tag "Create empty long-read evidence sentinels"

  input:
    path(empty_default_gtf)
    path(empty_alt_gtf)

  output:
    path "merged_transcriptomes.minimap2.long_reads.default_args.gtf", emit: default_args_gtf
    path "merged_transcriptomes.minimap2.long_reads.alt_args.gtf", emit: alt_args_gtf
    path "versions.yml", emit: versions

  script:
    """
    set -euo pipefail
    cp "${empty_default_gtf}" merged_transcriptomes.minimap2.long_reads.default_args.gtf
    cp "${empty_alt_gtf}" merged_transcriptomes.minimap2.long_reads.alt_args.gtf
    printf '"%s":\n  container: "not_recorded"\n' "${task.process}" > versions.yml
    """

  stub:
    """
    set -euo pipefail
    : > merged_transcriptomes.minimap2.long_reads.default_args.gtf
    : > merged_transcriptomes.minimap2.long_reads.alt_args.gtf
    printf '"%s":\n  container: "not_recorded"\n' "${task.process}" > versions.yml
    """
}
