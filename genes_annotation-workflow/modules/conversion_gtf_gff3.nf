// 5. Psiclass generates a GTF file, so we need to convert it in GFF3 format
process conversion_gtf_gff3 {

  tag "Convert GTF to GFF3 for transcriptome assembly"
  container 'avelt/gffread_gffutils:latest'
  containerOptions "--volume $params.outdir/evidence_data/transcriptomes:/transcriptomes --volume ${projectDir}/scripts/:/scripts --volume $genome_path:/genome_path"
  publishDir "$params.outdir/evidence_data/transcriptomes/rnaseq_stranded", pattern: "*RNAseq_stranded*"
  publishDir "$params.outdir/evidence_data/transcriptomes/rnaseq_unstranded", pattern: "*RNAseq_unstranded*"
  cpus 4

  input:
    val(genome_path)
    val(genome)
    path(assembly_transcriptome_stranded)
    path(assembly_transcriptome_unstranded)

  output:
    path("transcriptome_RNAseq_stranded.gff3"), emit : stranded_gff3
    path("transcriptome_RNAseq_stranded.gff3.transcripts.fasta"), emit : stranded_fasta
    path("transcriptome_RNAseq_unstranded.gff3"), emit : unstranded_gff3
    path("transcriptome_RNAseq_unstranded.gff3.transcripts.fasta"), emit : unstranded_fasta

  script:
    """
    /scripts/psiclass_gff3_formatting.sh -d /scripts -g /genome_path/$genome -i /transcriptomes/rnaseq_stranded/RNAseq_stranded_vote.gtf -o transcriptome_RNAseq_stranded.gff3 -f transcriptome_RNAseq_stranded.gff3.transcripts.fasta -s RNAseq_stranded
    /scripts/psiclass_gff3_formatting.sh -d /scripts -g /genome_path/$genome -i /transcriptomes/rnaseq_unstranded/RNAseq_unstranded_vote.gtf -o transcriptome_RNAseq_unstranded.gff3 -f transcriptome_RNAseq_unstranded.gff3.transcripts.fasta -s RNAseq_unstranded
    """
}
