process glimmerhmm_prediction {

  tag "Annotations prediction using GlimmerHMM"
  container 'avelt/glimmerhmm_gffutils:latest'
  containerOptions "--volume $params.outdir/GlimmerHMM/:/glimmerhmm"
  publishDir "$params.outdir/GlimmerHMM/"
  cpus 4

  input:
    val(genome_path)
    val(genome)
    val(stranded_gff3)

  output:
    path("training")

  script:
    """
    glimmerhmm {params.chromosome_file} {input.training_dir} -o glimmerhmm_prediction_${chromosome}.gff -n 1 -g
    """
}
