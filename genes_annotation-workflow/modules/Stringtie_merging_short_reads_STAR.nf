// Transcriptomes merging with StringTie
process Stringtie_merging_short_reads_STAR {

  tag "STAR/StringTie merging - short reads"
  container 'avelt/stringtie:latest'
  containerOptions "--volume ${projectDir}/scripts/:/scripts --volume ${projectDir}/work:/work --volume $params.outdir/evidence_data/transcriptomes/StringTie/short_reads/STAR/:/StringTie_short_reads_STAR"
  publishDir "$params.outdir/evidence_data/transcriptomes/StringTie/short_reads/STAR/"
  cpus 4

  input:
    val(concat_star_stringtie_annot)

  output:
    file("merged_transcriptomes.gtf")

  script:
    """
    gtf=\$(/scripts/retrieve_path_gtf.sh /StringTie_short_reads_STAR)
    /scripts/Stringtie_merging.sh -o merged_transcriptomes.gtf -g \${gtf}
    """
}
