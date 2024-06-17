process exonerate_mapping {

  tag "Exonerate alignment on the following protein sequences $organism"
  container 'avelt/exonerate_bedtools_samtools:latest'
  containerOptions "--volume $params.outdir/evidence_data/protein/:/exonerate --volume $params.outdir/evidence_data/protein/fasta_splitted:/fasta_splitted --volume ${projectDir}/scripts/:/scripts --volume $genome_path:/genome_path"
  publishDir "$params.outdir/evidence_data/protein_alignments_split/"
  cpus 4
  echo true

  input:
    val(genome_path)
    val(genome)
    tuple val(organism), path(psl_file)

  output:
    tuple val(organism), file("*.gff")

  script:
    """
    protein="\$(basename $psl_file | sed 's/.fasta.psl//')"
    /scripts/exonerate.sh -g /genome_path/$genome -a $psl_file -q /fasta_splitted/$organism/\${protein}.fasta -o \${protein} -d /scripts -e /exonerate
    """
}
