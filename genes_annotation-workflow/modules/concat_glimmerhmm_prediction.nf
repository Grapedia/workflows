process concat_glimmerhmm_prediction {

  tag "Concat all GlimmerHMM predictions"
  container 'avelt/glimmerhmm_gffutils:latest'
  containerOptions "--volume ${projectDir}/scripts:/scripts --volume ${projectDir}/work:/work --volume $params.outdir/GlimmerHMM/:/glimmerhmm"
  publishDir "$params.outdir/GlimmerHMM/"
  cpus 4

  input:
    val(glimmerhmm_pred)

  output:
    file("glimmerhmm_predictions.gff")

  script:
    """
    if ls /glimmerhmm/glimmerhmm_predictions.gff 1> /dev/null 2>&1
    then
      rm /glimmerhmm/glimmerhmm_predictions.gff
    fi
    gff=\$(/scripts/retrieve_path_gff.sh /glimmerhmm)
    cat \${gff} > glimmerhmm_predictions.gff
    """
}
