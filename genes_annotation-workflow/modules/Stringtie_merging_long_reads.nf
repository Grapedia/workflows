// Transcriptomes merging with StringTie
process Stringtie_merging_long_reads {

  tag "StringTie merging - long reads"
  container 'avelt/stringtie:latest'
  containerOptions "--volume ${projectDir}/scripts/:/scripts --volume ${projectDir}/work:/work --volume $params.outdir/evidence_data/transcriptomes/StringTie/long_reads/:/StringTie_long_reads"
  publishDir "$params.outdir/evidence_data/transcriptomes/StringTie/long_reads"
  cpus 4

  input:
    val(concat_minimap2_stringtie_annot)

  output:
    file("merged_transcriptomes.gtf")

  when:
  has_long_reads

  script:
    """
    gtf=\$(/scripts/retrieve_path_gtf.sh /StringTie_long_reads)
    /scripts/Stringtie_merging.sh -o merged_transcriptomes.gtf -g \${gtf}
    """
}
