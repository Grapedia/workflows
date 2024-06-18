// 4. Transcriptome assembly
process assembly_transcriptome_unstranded {

  tag "psiclass transcriptome assembly - unstranded"
  container 'avelt/psiclass_samtools:latest'
  containerOptions "--volume ${projectDir}/scripts:/scripts --volume ${projectDir}/work:/work --volume $params.outdir/evidence_data/RNAseq_unstranded/alignments/new_assembly:/alignments"
  publishDir "$params.outdir/evidence_data/transcriptomes/RNAseq_unstranded"
  cpus 4

  input:
    tuple val(sample_ID)

  output:
    file("RNAseq_unstranded_vote.gtf")

  shell:
    """
    bam=\$(/scripts/retrieve_path_bam.sh /alignments)
    /PsiCLASS-1.0.2/psiclass -p ${task.cpus} -b \${bam} -o RNAseq_unstranded
    """
}
