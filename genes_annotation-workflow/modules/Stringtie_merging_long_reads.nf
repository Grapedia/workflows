// Transcriptomes merging with StringTie
process Stringtie_merging_long_reads {

  tag "StringTie merging - long reads"
  container 'avelt/stringtie:latest'
  containerOptions "--volume ${projectDir}/scripts/:/scripts --volume ${projectDir}/work:/work --volume $params.outdir/evidence_data/transcriptomes/StringTie/long_reads/:/StringTie_long_reads"
  publishDir "$projectDir/FINAL_OUTPUT/transcriptomes/StringTie/long_reads"
  cpus 4

  input:
    val(concat_minimap2_stringtie_annot)

  output:
    file("*.gtf")

  when:
  has_long_reads

  script:
    """
    gtf_default=\$(/scripts/retrieve_path_gtf.sh /StringTie_long_reads default)
    gtf_alt=\$(/scripts/retrieve_path_gtf.sh /StringTie_long_reads alt)
    /scripts/Stringtie_merging.sh -o merged_transcriptomes.minimap2.long_reads.default_args.gtf -g \${gtf_default}
    /scripts/Stringtie_merging.sh -o merged_transcriptomes.minimap2.long_reads.alt_args.gtf -g \${gtf_alt}
    """
}
