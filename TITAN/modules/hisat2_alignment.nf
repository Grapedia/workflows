// 2. Aligning RNAseq data on reference genome with HISAT2
process hisat2_alignment {

  tag "HISAT2 on ${sample_ID} (${strand_type})"
  container 'avelt/hisat2:latest'
  containerOptions "--volume ${params.output_dir}/intermediate_files/evidence_data/hisat2_databases/:/hisat2_databases"

  publishDir { 
    if (strand_type == "unstranded") {
      return "${params.output_dir}/intermediate_files/evidence_data/RNAseq_alignments/HISAT2/unstranded"
    } else if (strand_type in ["stranded_forward", "stranded_reverse"]) {
      return "${params.output_dir}/intermediate_files/evidence_data/RNAseq_alignments/HISAT2/stranded"
    }
  }

  cpus 4

  input:
    val(hisat2_databases)
    tuple val(sample_ID), val(library_layout), path(reads), val(strand_type)
    val(genome)

  output:
    tuple val(sample_ID), path("${sample_ID}_Aligned.sort.bam"), val(strand_type), emit: samples_aligned

  script:
    """
    DATE=\$(date "+%Y-%m-%d %H:%M:%S")
    echo "[\$DATE] Running HISAT2 alignment on $sample_ID"
    if [[ $library_layout == "paired" ]]
    then
      if [[ "${strand_type}" == "unstranded" ]]; then
        CMD="/hisat2-2.2.1/hisat2 -p ${task.cpus} -x /hisat2_databases/$genome -1 ${reads[0]} -2 ${reads[1]} | samtools sort -o ${sample_ID}_Aligned.sort.bam -"
        echo "[\$DATE] Executing: \$CMD"
        /hisat2-2.2.1/hisat2 -p ${task.cpus} -x /hisat2_databases/$genome -1 ${reads[0]} -2 ${reads[1]} | samtools sort -o ${sample_ID}_Aligned.sort.bam -
      elif [[ "${strand_type}" == "stranded_forward" ]]; then
        CMD="/hisat2-2.2.1/hisat2 -p ${task.cpus} -x /hisat2_databases/$genome --rna-strandness FR -1 ${reads[0]} -2 ${reads[1]} | samtools sort -o ${sample_ID}_Aligned.sort.bam -"
        echo "[\$DATE] Executing: \$CMD"
        /hisat2-2.2.1/hisat2 -p ${task.cpus} -x /hisat2_databases/$genome --rna-strandness FR -1 ${reads[0]} -2 ${reads[1]} | samtools sort -o ${sample_ID}_Aligned.sort.bam -
      elif [[ "${strand_type}" == "stranded_reverse" ]]; then
        CMD="/hisat2-2.2.1/hisat2 -p ${task.cpus} -x /hisat2_databases/$genome --rna-strandness RF -1 ${reads[0]} -2 ${reads[1]} | samtools sort -o ${sample_ID}_Aligned.sort.bam -"
        echo "[\$DATE] Executing: \$CMD"
        /hisat2-2.2.1/hisat2 -p ${task.cpus} -x /hisat2_databases/$genome --rna-strandness RF -1 ${reads[0]} -2 ${reads[1]} | samtools sort -o ${sample_ID}_Aligned.sort.bam -
      fi
    elif [[ $library_layout == "single" ]]
    then
      if [[ "${strand_type}" == "unstranded" ]]; then
        CMD="/hisat2-2.2.1/hisat2 -p ${task.cpus} -x /hisat2_databases/$genome -U ${reads} | samtools sort -o ${sample_ID}_Aligned.sort.bam -"
        echo "[\$DATE] Executing: \$CMD"
        /hisat2-2.2.1/hisat2 -p ${task.cpus} -x /hisat2_databases/$genome -U ${reads} | samtools sort -o ${sample_ID}_Aligned.sort.bam -
      elif [[ "${strand_type}" == "stranded_forward" ]]; then
        CMD="/hisat2-2.2.1/hisat2 -p ${task.cpus} -x /hisat2_databases/$genome --rna-strandness FR -U ${reads} | samtools sort -o ${sample_ID}_Aligned.sort.bam -"
        echo "[\$DATE] Executing: \$CMD"
        /hisat2-2.2.1/hisat2 -p ${task.cpus} -x /hisat2_databases/$genome --rna-strandness FR -U ${reads} | samtools sort -o ${sample_ID}_Aligned.sort.bam -
      elif [[ "${strand_type}" == "stranded_reverse" ]]; then
        CMD="/hisat2-2.2.1/hisat2 -p ${task.cpus} -x /hisat2_databases/$genome --rna-strandness RF -U ${reads} | samtools sort -o ${sample_ID}_Aligned.sort.bam -"
        echo "[\$DATE] Executing: \$CMD"
        /hisat2-2.2.1/hisat2 -p ${task.cpus} -x /hisat2_databases/$genome --rna-strandness RF -U ${reads} | samtools sort -o ${sample_ID}_Aligned.sort.bam -
      fi
    fi
    """
}
