// Transcriptomes merging with StringTie
process Stringtie_merging_short_reads_STAR {

  tag "STAR/StringTie merging - short reads"
  container 'avelt/stringtie:latest'
  containerOptions "--volume ${projectDir}/scripts/:/scripts --volume ${projectDir}/work:/work --volume $params.outdir/evidence_data/transcriptomes/StringTie/short_reads/STAR/stranded/:/StringTie_short_reads_STAR_stranded --volume $params.outdir/evidence_data/transcriptomes/StringTie/short_reads/STAR/unstranded/:/StringTie_short_reads_STAR_unstranded"
  publishDir "$projectDir/FINAL_OUTPUT/transcriptomes/StringTie/short_reads/STAR/"
  cpus 4

  input:
    val(concat_star_stringtie_annot)

  output:
    file("*.gtf")

  script:
    """
    gtf_stranded_default=\$(/scripts/retrieve_path_gtf.sh /StringTie_short_reads_STAR_stranded default)
    gtf_unstranded_default=\$(/scripts/retrieve_path_gtf.sh /StringTie_short_reads_STAR_unstranded default)
    gtf_stranded_alt=\$(/scripts/retrieve_path_gtf.sh /StringTie_short_reads_STAR_stranded alt)
    gtf_unstranded_alt=\$(/scripts/retrieve_path_gtf.sh /StringTie_short_reads_STAR_unstranded alt)
    /scripts/Stringtie_merging.sh -o merged_transcriptomes.STAR.short_reads.default_args.stranded.gtf -g \${gtf_stranded_default}
    /scripts/Stringtie_merging.sh -o merged_transcriptomes.STAR.short_reads.default_args.unstranded.gtf -g \${gtf_unstranded_default}
    /scripts/Stringtie_merging.sh -o merged_transcriptomes.STAR.short_reads.alt_args.stranded.gtf -g \${gtf_stranded_alt}
    /scripts/Stringtie_merging.sh -o merged_transcriptomes.STAR.short_reads.alt_args.unstranded.gtf -g \${gtf_unstranded_alt}
    """
}
