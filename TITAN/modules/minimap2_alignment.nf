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

    CMD="minimap2 -t ${task.cpus} -ax splice:hq -uf ${minimap2_genome_indices} \${reads} | samtools view -b - | samtools sort -o ${sample_ID}_Aligned.sorted.bam -"
    echo "[\$DATE] Executing: \$CMD"
    minimap2 -t ${task.cpus} -ax splice:hq -uf ${minimap2_genome_indices} "\${reads}" \\
      | samtools view -b - \\
      | samtools sort -o ${sample_ID}_Aligned.sorted.bam -
    printf '"%s":\n  container: "not_recorded"\n' "${task.process}" > versions.yml
    """

  stub:
    """
    printf "BAM\\n" > ${sample_ID}_Aligned.sorted.bam
    printf '"%s":\n  container: "not_recorded"\n' "${task.process}" > versions.yml
    """
}
