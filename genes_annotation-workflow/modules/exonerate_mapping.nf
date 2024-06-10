process exonerate_mapping {

  tag "Exonerate alignment on the following protein sequences $organism"
  container 'quay.io/biocontainers/exonerate:2.4.0--h09da616_8'
  containerOptions "--volume $params.outdir/evidence_data/protein/fasta_splitted:/fasta_splitted --volume ${projectDir}/scripts/:/scripts --volume $genome_path:/genome_path"
  publishDir "$params.outdir/evidence_data/protein_alignments_split/"
  cpus 4

  input:
    val(genome_path)
    val(genome)
    tuple val(organism), path(psl_file)

  output:
    tuple val(organism), file("*.gff")

  script:
    protein = ${psl_file.baseName}.replaceFirst(/.psl/, "")
    """
    /scripts/exonerate.sh -g /genome_path/$genome -a $psl_file -q /fasta_splitted/$organism/$protein -o ${protein}.gff
    """
}
