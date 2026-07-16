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
    tuple val(sample_ID), val(library_layout), path(reads), val(strand_type)
    val(star_memory_per_job)

  output:
    tuple val(sample_ID), path("${sample_ID}_Aligned.sortedByCoord.out.bam"), val(strand_type), emit: samples_aligned
    path "versions.yml", emit: versions



  script:
    """
    set -euo pipefail
    DATE=\$(date "+%Y-%m-%d %H:%M:%S")
    echo "[\$DATE] Running STAR alignment on $sample_ID"

    if [[ $library_layout == "paired" ]]
    then

      CMD="STAR --readFilesCommand zcat --genomeDir ${star_database} --runThreadN ${task.cpus} --readFilesIn ${reads[0]} ${reads[1]} --outFileNamePrefix ${sample_ID}_ --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard --outSAMstrandField intronMotif --limitBAMsortRAM ${star_memory_per_job}"
      echo "[\$DATE] Executing: \$CMD"

      STAR --readFilesCommand zcat --genomeDir ${star_database} --runThreadN ${task.cpus} --readFilesIn ${reads[0]} ${reads[1]} --outFileNamePrefix ${sample_ID}_ --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard --outSAMstrandField intronMotif --limitBAMsortRAM ${star_memory_per_job}
    elif [[ $library_layout == "single" ]]
    then

      CMD="STAR --readFilesCommand zcat --genomeDir ${star_database} --runThreadN ${task.cpus} --readFilesIn ${reads} --outFileNamePrefix ${sample_ID}_ --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard --outSAMstrandField intronMotif --limitBAMsortRAM ${star_memory_per_job}"
      echo "[\$DATE] Executing: \$CMD"

      STAR --readFilesCommand zcat --genomeDir ${star_database} --runThreadN ${task.cpus} --readFilesIn ${reads} --outFileNamePrefix ${sample_ID}_ --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard --outSAMstrandField intronMotif --limitBAMsortRAM ${star_memory_per_job}
    fi
    printf '"%s":\n  container: "not_recorded"\n' "${task.process}" > versions.yml
    """

  stub:
    """
    printf "BAM\\n" > ${sample_ID}_Aligned.sortedByCoord.out.bam
    printf '"%s":\n  container: "not_recorded"\n' "${task.process}" > versions.yml
    """
}
