// 4. Transcriptome assembly with PsiCLASS
process assembly_transcriptome_star_psiclass {

  tag "STAR/PsiCLASS - short reads"
  container 'avelt/psiclass_samtools:latest'
  containerOptions "--volume ${projectDir}/scripts:/scripts --volume ${projectDir}/work:/work --volume $params.outdir/evidence_data/RNAseq_alignments/STAR:/alignments"
  publishDir "$params.outdir/evidence_data/transcriptomes/PsiCLASS"
  cpus 4

  input:
    val(sample_ID)

  output:
    file("RNAseq_vote.gtf")

  script:
    """
    bam=\$(/scripts/retrieve_path_bam.sh /alignments)
    /PsiCLASS-1.0.2/psiclass -p ${task.cpus} -b \${bam} -o RNAseq
    """
}
