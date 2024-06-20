process glimmerhmm_prediction {

  tag "Annotations prediction using GlimmerHMM"
  container 'avelt/glimmerhmm_gffutils:latest'
  containerOptions "--volume ${projectDir}/work:/work"
  publishDir "$params.outdir/GlimmerHMM/"
  cpus 4

  input:
    val(chr_fasta_file)
    val(glimmerhmm_training)

  output:
    file("*.gff")

  script:
    """
    chromosome=\$(basename $chr_fasta_file | sed "s/.fa//")
    chr_fasta_file_ok=\$(echo $chr_fasta_file | sed "s/.*work//")
    glimmerhmm_training_ok=\$(echo $glimmerhmm_training | sed "s/.*work//")
    glimmerhmm_linux_x86_64 /work/\$chr_fasta_file_ok /work/\$glimmerhmm_training_ok -o glimmerhmm_prediction_\$chromosome.gff -n 1 -g
    """
}
