// 5. Psiclass generates a GTF file, so we need to convert it in GFF3 format
process conversion_gtf_gff3 {

  tag "Convert GTF to GFF3 for transcriptome assembly"
  container 'quay.io/biocontainers/gffread:0.12.7--hdcf5f25_4'
  containerOptions "--volume $params.outdir/evidence_data/transcriptomes:/transcriptomes --volume ${projectDir}/scripts/:/scripts --volume $genome_path:/genome_path"
  publishDir "$params.outdir/evidence_data/transcriptomes/RNAseq_stranded", pattern: "*RNAseq_stranded*"
  publishDir "$params.outdir/evidence_data/transcriptomes/RNAseq_unstranded", pattern: "*RNAseq_unstranded*"
  cpus 4

  input:
    val(genome_path)
    val(genome)
    path(assembly_transcriptome_stranded)
    path(assembly_transcriptome_unstranded)

  output:
    file("transcriptome_RNAseq_stranded.gff3")
    file("transcriptome_RNAseq_stranded.gff3.transcripts.fasta")
    file("transcriptome_RNAseq_unstranded.gff3")
    file("transcriptome_RNAseq_unstranded.gff3.transcripts.fasta")

  script:
    """
    /scripts/psiclass_gff3_formatting.sh -d /scripts -g /genome_path/$genome -i /transcriptomes/RNAseq_stranded/RNAseq_stranded_vote.gtf -o transcriptome_RNAseq_stranded.gff3 -f transcriptome_RNAseq_stranded.gff3.transcripts.fasta -s RNAseq_stranded
    /scripts/psiclass_gff3_formatting.sh -d /scripts -g /genome_path/$genome -i /transcriptomes/RNAseq_unstranded/RNAseq_unstranded_vote.gtf -o transcriptome_RNAseq_unstranded.gff3 -f transcriptome_RNAseq_unstranded.gff3.transcripts.fasta -s RNAseq_unstranded
    """
}
