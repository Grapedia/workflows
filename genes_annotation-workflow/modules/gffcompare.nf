// 2. gffcompare
process gffcompare {

  tag "gffcompare on ${merged_gtf}"
  container 'avelt/gffcompare:latest'
  containerOptions "--volume ${projectDir}/work:/work --volume ${projectDir}/scripts:/scripts"
  publishDir "$params.outdir/evidence_data/gffcompare_HISAT2/"
  cpus 4

  input:
    val(merged_gtf)

  output:
    file("gffcmp*")

  script:
    """
    /scripts/retrieve_path_transcriptome_gffcompare.sh $merged_gtf
    /gffcompare-0.12.6/gffcompare -i gtf_list.txt
    """
}
