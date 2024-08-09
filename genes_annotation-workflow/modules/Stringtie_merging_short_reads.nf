// Transcriptomes merging with StringTie
process Stringtie_merging_short_reads {

  tag "StringTie merging - short reads"
  container 'avelt/stringtie:latest'
  containerOptions "--volume ${projectDir}/scripts/:/scripts --volume ${projectDir}/work:/work --volume $params.outdir/evidence_data/transcriptomes/StringTie/short_reads/:/StringTie_short_reads"
  publishDir "$params.outdir/evidence_data/transcriptomes/StringTie/short_reads"
  cpus 4

  input:
    val(concat_star_stringtie_annot)

  output:
    file("merged_transcriptomes.gtf")

  script:
    """
    gtf=\$(/scripts/retrieve_path_gtf.sh /StringTie_short_reads)
    /scripts/Stringtie_merging.sh -o merged_transcriptomes.gtf -g \${gtf}
    """
}
