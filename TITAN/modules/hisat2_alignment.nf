// 2. Aligning RNAseq data on reference genome with HISAT2
process hisat2_alignment {
  label 'process_alignment'

  tag "HISAT2 on ${sample_ID} (${strand_type})"
  container params.container_hisat2

  publishDir {
    if (strand_type == "unstranded") {
      return "${params.output_dir}/intermediate_files/evidence_data/RNAseq_alignments/HISAT2/unstranded"
    } else if (strand_type in ["stranded_forward", "stranded_reverse"]) {
      return "${params.output_dir}/intermediate_files/evidence_data/RNAseq_alignments/HISAT2/stranded"
    }
  }, mode: "copy", enabled: params.publish_intermediates
  input:
    path(hisat2_databases)
    tuple val(sample_ID), val(library_layout), path(read_1), path(read_2), val(strand_type)
    val(genome)

  output:
    tuple val(sample_ID), path("${sample_ID}_Aligned.sort.bam"), val(strand_type), emit: samples_aligned
    path "versions.yml", emit: versions



  script:
    """
    set -euo pipefail
    DATE=\$(date "+%Y-%m-%d %H:%M:%S")
    echo "[\$DATE] Running HISAT2 alignment on $sample_ID"
    HISAT2_INDEX="${hisat2_databases}/${genome}"

    if [[ $library_layout == "paired" ]]
    then
      if [[ "${strand_type}" == "unstranded" ]]; then
        CMD="/hisat2-2.2.1/hisat2 -p ${task.cpus} -x \${HISAT2_INDEX} -1 ${read_1} -2 ${read_2} | samtools sort -o ${sample_ID}_Aligned.sort.bam -"
        echo "[\$DATE] Executing: \$CMD"
        /hisat2-2.2.1/hisat2 -p ${task.cpus} -x "\${HISAT2_INDEX}" -1 ${read_1} -2 ${read_2} | samtools sort -o ${sample_ID}_Aligned.sort.bam -
      elif [[ "${strand_type}" == "stranded_forward" ]]; then
        CMD="/hisat2-2.2.1/hisat2 -p ${task.cpus} -x \${HISAT2_INDEX} --rna-strandness FR -1 ${read_1} -2 ${read_2} | samtools sort -o ${sample_ID}_Aligned.sort.bam -"
        echo "[\$DATE] Executing: \$CMD"
        /hisat2-2.2.1/hisat2 -p ${task.cpus} -x "\${HISAT2_INDEX}" --rna-strandness FR -1 ${read_1} -2 ${read_2} | samtools sort -o ${sample_ID}_Aligned.sort.bam -
      elif [[ "${strand_type}" == "stranded_reverse" ]]; then
        CMD="/hisat2-2.2.1/hisat2 -p ${task.cpus} -x \${HISAT2_INDEX} --rna-strandness RF -1 ${read_1} -2 ${read_2} | samtools sort -o ${sample_ID}_Aligned.sort.bam -"
        echo "[\$DATE] Executing: \$CMD"
        /hisat2-2.2.1/hisat2 -p ${task.cpus} -x "\${HISAT2_INDEX}" --rna-strandness RF -1 ${read_1} -2 ${read_2} | samtools sort -o ${sample_ID}_Aligned.sort.bam -
      fi
    elif [[ $library_layout == "single" ]]
    then
      if [[ "${strand_type}" == "unstranded" ]]; then
        CMD="/hisat2-2.2.1/hisat2 -p ${task.cpus} -x \${HISAT2_INDEX} -U ${read_1} | samtools sort -o ${sample_ID}_Aligned.sort.bam -"
        echo "[\$DATE] Executing: \$CMD"
        /hisat2-2.2.1/hisat2 -p ${task.cpus} -x "\${HISAT2_INDEX}" -U ${read_1} | samtools sort -o ${sample_ID}_Aligned.sort.bam -
      elif [[ "${strand_type}" == "stranded_forward" ]]; then
        CMD="/hisat2-2.2.1/hisat2 -p ${task.cpus} -x \${HISAT2_INDEX} --rna-strandness FR -U ${read_1} | samtools sort -o ${sample_ID}_Aligned.sort.bam -"
        echo "[\$DATE] Executing: \$CMD"
        /hisat2-2.2.1/hisat2 -p ${task.cpus} -x "\${HISAT2_INDEX}" --rna-strandness FR -U ${read_1} | samtools sort -o ${sample_ID}_Aligned.sort.bam -
      elif [[ "${strand_type}" == "stranded_reverse" ]]; then
        CMD="/hisat2-2.2.1/hisat2 -p ${task.cpus} -x \${HISAT2_INDEX} --rna-strandness RF -U ${read_1} | samtools sort -o ${sample_ID}_Aligned.sort.bam -"
        echo "[\$DATE] Executing: \$CMD"
        /hisat2-2.2.1/hisat2 -p ${task.cpus} -x "\${HISAT2_INDEX}" --rna-strandness RF -U ${read_1} | samtools sort -o ${sample_ID}_Aligned.sort.bam -
      fi
    fi
    printf '"%s":\n  container: "not_recorded"\n' "${task.process}" > versions.yml
    """

  stub:
    """
    printf "BAM\\n" > ${sample_ID}_Aligned.sort.bam
    printf '"%s":\n  container: "not_recorded"\n' "${task.process}" > versions.yml
    """
}
