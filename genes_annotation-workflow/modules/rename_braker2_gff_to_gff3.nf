process rename_braker2_gff_to_gff3 {

  tag "Executing BRAKER2 output gff3 renaming"
  container 'avelt/braker2_prothint_genemark:latest'
  publishDir "$params.outdir/BRAKER2"
  cpus 4

  input:
    path(braker2_prediction)

  output:
    path("braker2_predictions.gff3")

  script:
    """
    mv $braker2_prediction braker2_predictions.gff3
    """
}
