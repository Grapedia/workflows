// 5. Psiclass generates a GTF file, so we need to convert it in GFF3 format
process conversion_gtf_gff3_star_psiclass {

  tag "Convert GTF to GFF3 for STAR/PsiCLASS transcriptome assembly"
  container 'avelt/gffread_gffutils:latest'
  containerOptions "--volume $params.outdir/evidence_data/transcriptomes/PsiCLASS:/transcriptomes --volume ${projectDir}/scripts/:/scripts --volume $genome_path:/genome_path"
  publishDir "$params.outdir/evidence_data/transcriptomes/PsiCLASS"
  cpus 4

  input:
    val(genome_path)
    val(genome)
    path(assembly_transcriptome_star_psiclass)

  output:
    path("transcriptome_RNAseq_star_psiclass.gff3"), emit : gff3
    path("transcriptome_RNAseq_star_psiclass.gff3.transcripts.fasta"), emit : fasta

  script:
    """
    gtf=\$(/scripts/retrieve_path_gtf.sh /transcriptomes)
    /scripts/psiclass_gff3_formatting.sh -d /scripts -g /genome_path/$genome -i \${gtf} -o transcriptome_RNAseq_star_psiclass.gff3 -f transcriptome_RNAseq_star_psiclass.gff3.transcripts.fasta
    """
}
