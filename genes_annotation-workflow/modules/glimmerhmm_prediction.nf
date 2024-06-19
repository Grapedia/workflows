process glimmerhmm_prediction {

  tag "Annotations prediction using GlimmerHMM"
  container 'avelt/glimmerhmm_gffutils:latest'
  publishDir "$params.outdir/GlimmerHMM/"
  cpus 4

  input:
    val(chr_fasta_file)
    val(glimmerhmm_training)

  output:
    file("*.gff")

  script:
    """
    chromosome=\$(echo \$chr_fasta_file | basename)
    echo \$chromosome
    glimmerhmm $chr_fasta_file $glimmerhmm_training -o glimmerhmm_prediction_\$chromosome.gff -n 1 -g
    """
}
