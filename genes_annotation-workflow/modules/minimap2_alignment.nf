// 2. Aligning RNAseq data on reference genome with minimap2
process minimap2_alignment {

  tag "Minimap2 on ${sample_ID}"
  container 'avelt/minimap2_samtools:latest'
  containerOptions "--volume ${projectDir}/data/RNAseq_data:/RNAseq_data --volume ${projectDir}/work:/work --volume ${projectDir}/scripts:/scripts"
  publishDir "$params.outdir/evidence_data/RNAseq_alignments/minimap2"
  cpus 4

  input:
    val(minimap2_genome_indices)
    tuple val(sample_ID), val(SRA_or_FASTQ), val(library_layout)

  output:
    file("${sample_ID}_Aligned.sorted.bam")

  script:
    // PacBio Iso-seq/traditional cDNA parameters
    """
    minimap2_index=\$(/scripts/retrieve_path_minimap2_index.sh ${minimap2_genome_indices})
    minimap2 -t ${task.cpus} -ax splice:hq -uf \${minimap2_index} /RNAseq_data/${sample_ID}.fastq.gz > ${sample_ID}_Aligned.sam
    samtools view -b ${sample_ID}_Aligned.sam | samtools sort - > ${sample_ID}_Aligned.sorted.bam
    rm ${sample_ID}_Aligned.sam
    """
}
