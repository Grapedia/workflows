// 2. Aligning RNAseq data on reference genome with STAR
process star_alignment {
  label 'process_alignment'

  tag "STAR on ${sample_ID} (${strand_type})"
  container params.container_star
  publishDir {
    if (strand_type == "unstranded") {
      return "${params.output_dir}/intermediate_files/evidence_data/RNAseq_alignments/STAR/unstranded"
    } else if (strand_type in ["stranded_forward", "stranded_reverse"]) {
      return "${params.output_dir}/intermediate_files/evidence_data/RNAseq_alignments/STAR/stranded"
    }
  }, mode: "copy", enabled: params.publish_intermediates
  input:
    path(star_database)
    tuple val(sample_ID), val(library_layout), path(read_1), path(read_2), val(strand_type)

  output:
    // Downstream transcript assemblers consume the coordinate-sorted BAM directly;
    // no BAM index is required by the current STAR -> StringTie/PsiCLASS/BRAKER3 graph.
    tuple val(sample_ID), path("${sample_ID}_Aligned.sortedByCoord.out.bam"), val(strand_type), emit: samples_aligned
    path "versions.yml", emit: versions



  script:
    """
    set -euo pipefail
    DATE=\$(date "+%Y-%m-%d %H:%M:%S")
    echo "[\$DATE] Running STAR alignment on $sample_ID"
    LIMIT_BAM_SORT_RAM="${task.memory.toBytes()}"

    if [[ $library_layout == "paired" ]]; then
      CMD="STAR --readFilesCommand zcat --genomeDir ${star_database} --runThreadN ${task.cpus} --readFilesIn ${read_1} ${read_2} --outFileNamePrefix ${sample_ID}_ --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard --outSAMstrandField intronMotif --limitBAMsortRAM \${LIMIT_BAM_SORT_RAM}"
      echo "[\$DATE] Executing: \$CMD"

      STAR --readFilesCommand zcat --genomeDir ${star_database} --runThreadN ${task.cpus} --readFilesIn ${read_1} ${read_2} --outFileNamePrefix ${sample_ID}_ --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard --outSAMstrandField intronMotif --limitBAMsortRAM "\${LIMIT_BAM_SORT_RAM}"
    elif [[ $library_layout == "single" ]]; then
      CMD="STAR --readFilesCommand zcat --genomeDir ${star_database} --runThreadN ${task.cpus} --readFilesIn ${read_1} --outFileNamePrefix ${sample_ID}_ --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard --outSAMstrandField intronMotif --limitBAMsortRAM \${LIMIT_BAM_SORT_RAM}"
      echo "[\$DATE] Executing: \$CMD"

      STAR --readFilesCommand zcat --genomeDir ${star_database} --runThreadN ${task.cpus} --readFilesIn ${read_1} --outFileNamePrefix ${sample_ID}_ --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard --outSAMstrandField intronMotif --limitBAMsortRAM "\${LIMIT_BAM_SORT_RAM}"
    else
      echo "[\$DATE] ERROR: library_layout must be 'paired' or 'single', got '$library_layout'" >&2
      exit 1
    fi

    STAR --version 2>&1 | sed 's/^/  star: "/; s/\$/"/' | {
      printf '"%s":\\n  limitBAMsortRAM: "%s"\\n' "${task.process}" "\$LIMIT_BAM_SORT_RAM"
      cat
    } > versions.yml
    """

  stub:
    """
    set -euo pipefail

    printf "BAM\\n" > ${sample_ID}_Aligned.sortedByCoord.out.bam
    printf '"%s":\n  star: "stub"\n  limitBAMsortRAM: "stub"\n' "${task.process}" > versions.yml
    """
}
