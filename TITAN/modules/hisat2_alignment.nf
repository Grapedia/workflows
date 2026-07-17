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
    path("${sample_ID}_Aligned.sort.bam.bai"), emit: bam_index
    path "versions.yml", emit: versions



  script:
    """
    set -euo pipefail
    DATE=\$(date "+%Y-%m-%d %H:%M:%S")
    echo "[\$DATE] Running HISAT2 alignment on $sample_ID"
    HISAT2_INDEX="${hisat2_databases}/${genome}"
    BAM="${sample_ID}_Aligned.sort.bam"
    HISAT2_STRANDNESS=""

    if [[ "${library_layout}" == "paired" ]]; then
      case "${strand_type}" in
        unstranded) HISAT2_STRANDNESS="" ;;
        stranded_forward) HISAT2_STRANDNESS="--rna-strandness FR" ;;
        stranded_reverse) HISAT2_STRANDNESS="--rna-strandness RF" ;;
        *)
          echo "Unsupported strand_type for HISAT2 paired-end alignment: ${strand_type}" >&2
          exit 1
          ;;
      esac
      CMD="/hisat2-2.2.1/hisat2 -p ${task.cpus} -x \${HISAT2_INDEX} \${HISAT2_STRANDNESS} -1 ${read_1} -2 ${read_2} | samtools sort -o \${BAM} -"
      echo "[\$DATE] Executing: \$CMD"
      /hisat2-2.2.1/hisat2 -p ${task.cpus} -x "\${HISAT2_INDEX}" \${HISAT2_STRANDNESS} -1 ${read_1} -2 ${read_2} | samtools sort -o "\${BAM}" -
    elif [[ "${library_layout}" == "single" ]]; then
      case "${strand_type}" in
        unstranded) HISAT2_STRANDNESS="" ;;
        stranded_forward) HISAT2_STRANDNESS="--rna-strandness F" ;;
        stranded_reverse) HISAT2_STRANDNESS="--rna-strandness R" ;;
        *)
          echo "Unsupported strand_type for HISAT2 single-end alignment: ${strand_type}" >&2
          exit 1
          ;;
      esac
      CMD="/hisat2-2.2.1/hisat2 -p ${task.cpus} -x \${HISAT2_INDEX} \${HISAT2_STRANDNESS} -U ${read_1} | samtools sort -o \${BAM} -"
      echo "[\$DATE] Executing: \$CMD"
      /hisat2-2.2.1/hisat2 -p ${task.cpus} -x "\${HISAT2_INDEX}" \${HISAT2_STRANDNESS} -U ${read_1} | samtools sort -o "\${BAM}" -
    else
      echo "Unsupported library_layout for HISAT2 alignment: ${library_layout}" >&2
      exit 1
    fi

    samtools index -@ ${task.cpus} "\${BAM}"
    {
      printf '"%s":\n' "${task.process}"
      /hisat2-2.2.1/hisat2 --version 2>&1 | awk 'NR == 1 { printf "  hisat2: \\"%s\\"\\n", \$0 }'
      samtools --version | awk 'NR == 1 { printf "  samtools: \\"%s\\"\\n", \$0 }'
    } > versions.yml
    """

  stub:
    """
    set -euo pipefail
    printf "BAM\\n" > ${sample_ID}_Aligned.sort.bam
    printf "BAI\\n" > ${sample_ID}_Aligned.sort.bam.bai
    printf '"%s":\n  hisat2: "stub"\n  samtools: "stub"\n' "${task.process}" > versions.yml
    """
}
