// To parallelize GlimmerHMM prediction, we split the assembly FASTA file in separate file ...
// ... (one for each chromosome) to run the prediction separatly on each chromosome
process split_fasta {

  tag "Executing fastaexplode on the new assembly"
  containerOptions "--volume $genome_path:/genome_path"
  container 'avelt/exonerate_bedtools_samtools:latest'
  publishDir "$params.outdir/GlimmerHMM/fastaexplode/"
  cpus 4

  input:
    val(genome_path)
    val(genome)

  output:
    file("*.fa")

  script:
    """
    fastaexplode -f /genome_path/${genome}
    """
}
