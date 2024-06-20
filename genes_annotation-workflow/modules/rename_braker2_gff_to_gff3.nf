process rename_braker2_gff_to_gff3 {

  tag "Executing BRAKER2 output gff3 renaming"
  container 'avelt/braker2_prothint_genemark:latest'
  publishDir "$params.outdir/BRAKER2_RNAseq_stranded", pattern: "*RNAseq_stranded*"
  publishDir "$params.outdir/BRAKER2_RNAseq_unstranded", pattern: "*RNAseq_unstranded*"
  cpus 4

  input:
    path(braker2_prediction_stranded)
    path(braker2_prediction_unstranded)

  output:
    path("braker2_predictions_RNAseq_stranded.gff3"), emit : braker2_prediction_stranded
    path("braker2_predictions_RNAseq_unstranded.gff3"), emit : braker2_prediction_unstranded

  script:
    """
    mv $braker2_prediction_stranded braker2_predictions_RNAseq_stranded.gff3
    mv $braker2_prediction_unstranded braker2_predictions_RNAseq_unstranded.gff3
    """
}
