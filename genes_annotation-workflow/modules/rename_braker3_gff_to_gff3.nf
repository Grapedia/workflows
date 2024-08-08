process rename_braker3_gff_to_gff3 {

  tag "Executing BRAKER3 output gff3 renaming"
  container 'avelt/braker3:latest'
  publishDir "$params.outdir/BRAKER3"
  cpus 4

  input:
    path(braker3_prediction)

  output:
    path("braker3_prediction.gff3")

  script:
    """
    mv $braker3_prediction braker3_prediction.gff3
    """
}
