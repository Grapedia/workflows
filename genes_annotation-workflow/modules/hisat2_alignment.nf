// 2. Aligning RNAseq data on reference genome with HISAT2
process hisat2_alignment {

  tag "HISAT2 on ${sample_ID}"
  container 'avelt/hisat2:latest'
  containerOptions "--volume $params.outdir/evidence_data/hisat2_databases/:/hisat2_databases"
  publishDir "$params.outdir/evidence_data/RNAseq_alignments/HISAT2"
  cpus 4

  input:
    val(hisat2_databases)
    tuple val(sample_ID), val(library_layout), path(reads)
    val(genome)

  output:
    file("${sample_ID}_Aligned.sort.bam")

  script:
    """
    if [[ $library_layout == "paired" ]]
    then
      /hisat2-2.2.1/hisat2 -p ${task.cpus} -x /hisat2_databases/$genome -1 ${reads[0]} -2 ${reads[1]} | samtools sort -o ${sample_ID}_Aligned.sort.bam -
    elif [[ $library_layout == "single" ]]
    then
      /hisat2-2.2.1/hisat2 -p ${task.cpus} -x /hisat2_databases/$genome -U ${reads} | samtools sort -o ${sample_ID}_Aligned.sort.bam -
    fi
    """
}
