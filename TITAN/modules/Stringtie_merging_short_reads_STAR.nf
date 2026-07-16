// Transcriptomes merging with StringTie
process Stringtie_merging_short_reads_STAR {
  label 'process_merge'

  tag "STAR/StringTie merge: stranded and optional unstranded short-read GTFs"
  container params.container_stringtie
  publishDir "${params.output_dir}/tmp", mode: 'copy', enabled: params.publish_intermediates
  publishDir "${params.output_dir}", mode: 'copy', saveAs: { filename ->
    if (filename == 'merged_transcriptomes.STAR.short_reads.default_args.stranded.gtf') {
      return 'merged_star_stringtie_stranded_default.gtf'
    }
    if (filename == 'merged_transcriptomes.STAR.short_reads.alt_args.stranded.gtf') {
      return 'merged_star_stringtie_stranded_alt.gtf'
    }
    if (filename == 'merged_transcriptomes.STAR.short_reads.default_args.unstranded.gtf') {
      return 'merged_star_stringtie_unstranded_default.gtf'
    }
    if (filename == 'merged_transcriptomes.STAR.short_reads.alt_args.unstranded.gtf') {
      return 'merged_star_stringtie_unstranded_alt.gtf'
    }
    return null
  }

  input:
    path(stranded_default_gtfs)
    path(stranded_alt_gtfs)
    path(unstranded_default_gtfs)
    path(unstranded_alt_gtfs)

  output:
    path "merged_transcriptomes.STAR.short_reads.default_args.stranded.gtf", emit: default_args_stranded
    path "merged_transcriptomes.STAR.short_reads.alt_args.stranded.gtf", emit: alt_args_stranded
    path "merged_transcriptomes.STAR.short_reads.default_args.unstranded.gtf", optional: true, emit: default_args_unstranded
    path "merged_transcriptomes.STAR.short_reads.alt_args.unstranded.gtf", optional: true, emit: alt_args_unstranded
    path "versions.yml", emit: versions



  script:
    """
    set -euo pipefail
    DATE=\$(date "+%Y-%m-%d %H:%M:%S")
    echo "[\$DATE] Running StringTie merging on STAR/StringTie transcriptomes - separating stranded and unstranded samples."

    printf '%s\\n' ${stranded_default_gtfs} > stranded_default_gtfs.txt
    printf '%s\\n' ${stranded_alt_gtfs} > stranded_alt_gtfs.txt
    stringtie --merge -o merged_transcriptomes.STAR.short_reads.default_args.stranded.gtf stranded_default_gtfs.txt
    stringtie --merge -o merged_transcriptomes.STAR.short_reads.alt_args.stranded.gtf stranded_alt_gtfs.txt

    unstranded_default_files=( ${unstranded_default_gtfs} )
    unstranded_alt_files=( ${unstranded_alt_gtfs} )
    if [[ \${#unstranded_default_files[@]} -gt 0 && \${#unstranded_alt_files[@]} -gt 0 && -s "\${unstranded_default_files[0]}" && -s "\${unstranded_alt_files[0]}" ]]; then
      echo "[\$DATE] Running StringTie merging on STAR/StringTie transcriptomes - unstranded samples detected."
      printf '%s\\n' ${unstranded_default_gtfs} > unstranded_default_gtfs.txt
      printf '%s\\n' ${unstranded_alt_gtfs} > unstranded_alt_gtfs.txt
      stringtie --merge -o merged_transcriptomes.STAR.short_reads.default_args.unstranded.gtf unstranded_default_gtfs.txt
      stringtie --merge -o merged_transcriptomes.STAR.short_reads.alt_args.unstranded.gtf unstranded_alt_gtfs.txt
    else
      : > merged_transcriptomes.STAR.short_reads.default_args.unstranded.gtf
      : > merged_transcriptomes.STAR.short_reads.alt_args.unstranded.gtf
    fi
    printf '"%s":\n  container: "not_recorded"\n' "${task.process}" > versions.yml
    """

  stub:
    """
    printf "chr1\\tStringTie\\ttranscript\\t1\\t10\\t.\\t+\\t.\\tgene_id \\"star_stranded_gene\\"; transcript_id \\"star_stranded_tx\\";\\n" > merged_transcriptomes.STAR.short_reads.default_args.stranded.gtf
    printf "chr1\\tStringTie\\ttranscript\\t1\\t10\\t.\\t+\\t.\\tgene_id \\"star_stranded_alt_gene\\"; transcript_id \\"star_stranded_alt_tx\\";\\n" > merged_transcriptomes.STAR.short_reads.alt_args.stranded.gtf
    printf "chr1\\tStringTie\\ttranscript\\t1\\t10\\t.\\t+\\t.\\tgene_id \\"star_unstranded_gene\\"; transcript_id \\"star_unstranded_tx\\";\\n" > merged_transcriptomes.STAR.short_reads.default_args.unstranded.gtf
    printf "chr1\\tStringTie\\ttranscript\\t1\\t10\\t.\\t+\\t.\\tgene_id \\"star_unstranded_alt_gene\\"; transcript_id \\"star_unstranded_alt_tx\\";\\n" > merged_transcriptomes.STAR.short_reads.alt_args.unstranded.gtf
    printf '"%s":\n  container: "not_recorded"\n' "${task.process}" > versions.yml
    """
}
