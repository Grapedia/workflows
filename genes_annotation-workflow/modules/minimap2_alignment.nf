// 2. Aligning RNAseq data on reference genome with minimap2
process minimap2_alignment {

  tag "Minimap2 on ${sample_ID}"
  container 'avelt/minimap2_samtools:latest'
  containerOptions "--memory=50g --volume $genome_path:/genome_path --volume ${projectDir}/data/RNAseq_data:/RNAseq_data"
  publishDir "$params.outdir/evidence_data/RNAseq_alignments/minimap2"
  cpus 4

  input:
    val(genome_path)
    val(genome)
    tuple val(sample_ID), val(SRA_or_FASTQ), val(library_layout)

  output:
    file("${sample_ID}_Aligned.sorted.bam")

  script:
    // PacBio Iso-seq/traditional cDNA parameters
    """
    minimap2 -d /genome_path/${genome}.mmi /genome_path/$genome
    minimap2 -t ${task.cpus} -ax splice:hq -uf /genome_path/${genome}.mmi /RNAseq_data/${sample_ID}.fastq.gz > ${sample_ID}_Aligned.sam
    samtools view -b ${sample_ID}_Aligned.sam | samtools sort - > ${sample_ID}_Aligned.sorted.bam
    rm ${sample_ID}_Aligned.sam
    """
}
