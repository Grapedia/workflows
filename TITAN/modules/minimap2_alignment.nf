// 2. Aligning RNAseq data on reference genome with minimap2
process minimap2_alignment {
  label 'process_alignment'

  tag "Minimap2 on ${sample_ID}"
  container params.container_minimap2_samtools
  publishDir "${params.output_dir}/intermediate_files/evidence_data/RNAseq_alignments/minimap2", mode: "copy", enabled: params.publish_intermediates
  input:
    path(minimap2_genome_indices)
    tuple val(sample_ID), val(SRA_or_FASTQ), val(library_layout), val(read_format), path(reads_fastq), path(reads_fasta)

  output:
    tuple val(sample_ID), path("${sample_ID}_Aligned.sorted.bam"), emit: samples_aligned
    path "versions.yml", emit: versions

  script:
    """
    set -euo pipefail
    DATE=\$(date "+%Y-%m-%d %H:%M:%S")
    echo "[\$DATE] Running minimap2 alignment of $sample_ID"

    if [[ $SRA_or_FASTQ != "FASTQ" && $SRA_or_FASTQ != "SRA" && $SRA_or_FASTQ != "FASTA" ]]; then
      echo "[\$DATE] ERROR: $SRA_or_FASTQ is not equal to FASTQ, SRA or FASTA" >&2
      exit 1
    fi

    if [[ "${read_format}" == "fastq" ]]; then
      reads="${reads_fastq}"
    elif [[ "${read_format}" == "fasta" ]]; then
      reads="${reads_fasta}"
    else
      echo "[\$DATE] ERROR: read_format must be 'fastq' or 'fasta', got '${read_format}'" >&2
      exit 1
    fi

    # samtools sort defaults to 1 thread / 768M buffer; give it a fair share
    # of task.cpus so it doesn't stall the minimap2 pipe (see hisat2_alignment.nf).
    TOTAL_CPUS=${task.cpus}
    if [ "\$TOTAL_CPUS" -ge 8 ]; then
      SORT_CPUS=4
    elif [ "\$TOTAL_CPUS" -ge 4 ]; then
      SORT_CPUS=2
    else
      SORT_CPUS=1
    fi
    MINIMAP2_CPUS=\$(( TOTAL_CPUS - SORT_CPUS ))
    if [ "\$MINIMAP2_CPUS" -lt 1 ]; then
      MINIMAP2_CPUS=1
    fi
    SORT_MEM="2G"

    CMD="minimap2 -t \${MINIMAP2_CPUS} -ax splice:hq -uf ${minimap2_genome_indices} \${reads} | samtools sort -@ \${SORT_CPUS} -m \${SORT_MEM} -o ${sample_ID}_Aligned.sorted.bam -"
    echo "[\$DATE] Executing: \$CMD"
    minimap2 -t "\${MINIMAP2_CPUS}" -ax splice:hq -uf ${minimap2_genome_indices} "\${reads}" \\
      | samtools sort -@ "\${SORT_CPUS}" -m "\${SORT_MEM}" -o ${sample_ID}_Aligned.sorted.bam -
    {
      printf '"%s":\n' "${task.process}"
      minimap2 --version | awk '{ printf "  minimap2: \\"%s\\"\\n", \$0 }'
      samtools --version | awk 'NR == 1 { printf "  samtools: \\"%s\\"\\n", \$0 }'
    } > versions.yml
    """

  stub:
    """
    set -euo pipefail
    printf "BAM\\n" > ${sample_ID}_Aligned.sorted.bam
    printf '"%s":\n  minimap2: "stub"\n  samtools: "stub"\n' "${task.process}" > versions.yml
    """
}
