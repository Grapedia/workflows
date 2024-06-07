process pblat_protein_alignment {

  tag "pblat alignment of ${filename} proteins against ${genome} genome"
  container 'quay.io/biocontainers/pblat:2.5--h0e0aaa8_2'
  containerOptions "--volume $params.outdir/evidence_data/protein/fasta_splitted/${organism}:/fasta_splitted_${organism} --volume $genome_path:/genome_path"
  publishDir "$params.outdir/evidence_data/protein/alignments_pblat/${organism}"
  cpus 4

  input:
    val(genome_path)
    val(genome)
    tuple val(organism), path(filename)

  output:
    tuple val(organism), file("*.psl")

  script:
    """
    pblat -threads=${task.cpus} -t=dnax -q=prot -minIdentity=90 -mask=lower -noHead /genome_path/$genome $filename ${filename}.psl
    """
}
