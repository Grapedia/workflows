// 4. Transcriptome assembly
process assembly_transcriptome_unstranded {

  tag "psiclass transcriptome assembly - unstranded"
  container 'avelt/psiclass_samtools:latest'
  containerOptions "--volume ${projectDir}/work:/work --volume $params.outdir/evidence_data/RNAseq_unstranded/alignments/new_assembly:/alignments"
  publishDir "$params.outdir/evidence_data/transcriptomes/RNAseq_unstranded"
  cpus 4

  input:
    tuple val(sample_ID)

  output:
    file("RNAseq_unstranded_vote.gtf")

  script:
    """
    bam="\$(ll /alignments/*bam | sed 's/.*-> .*work\//\/work\//' | tr '\n' ',' | sed 's/,\$//')"
    /PsiCLASS-1.0.2/psiclass -p ${task.cpus} -b \${bam} -o RNAseq_unstranded
    """
}
