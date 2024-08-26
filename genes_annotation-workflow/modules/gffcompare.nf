// 2. gffcompare
process gffcompare {

  tag "gffcompare on ${merged_gtf}"
  container 'avelt/gffcompare:latest'
  containerOptions ""
  publishDir "$params.outdir/evidence_data/gffcompare_HISAT2/"
  cpus 4

  input:
    val(merged_gtf)

  output:
    file("gffcmp*")

  script:
    """
    /gffcompare-0.12.6/gffcompare -i $merged_gtf
    """
}
