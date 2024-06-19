process rename_maker_gff_to_gff3 {

  tag "Executing MAKER output gff conversion to gff3"
  container 'quay.io/biocontainers/maker:3.01.03--pl526hb8757ab_0'
  publishDir "$params.outdir/MAKER"
  cpus 4

  input:
    val(maker_gff)

  output:
    file("maker_snap_predictions.gff3")

  script:
    """
    mv $maker_gff maker_snap_predictions.gff3
    """
}
