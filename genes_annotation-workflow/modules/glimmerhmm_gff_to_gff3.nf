process glimmerhmm_gff_to_gff3 {

  tag "Executing GlimmerHMM output gff conversion to gff3"
  container 'quay.io/biocontainers/evidencemodeler:1.1.1--0'
  publishDir "$params.outdir/GlimmerHMM/"
  cpus 4

  input:
    val(glimmerhmm_predictions)

  output:
    file("glimmerhmm_predictions.gff3")

  script:
    """
    /usr/local/opt/evidencemodeler-1.1.1/EvmUtils/misc/glimmerHMM_to_GFF3.pl $glimmerhmm_predictions > glimmerhmm_predictions.gff3
    """
}
