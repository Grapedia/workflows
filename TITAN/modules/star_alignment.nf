// 2. Aligning RNAseq data on reference genome with STAR
process star_alignment {

  tag "STAR on ${sample_ID} (${strand_type})"
  container 'quay.io/biocontainers/star:2.7.11b--h43eeafb_2'
  containerOptions "--memory=50g --volume $params.outdir/evidence_data/star_databases/:/star_databases"
  publishDir { 
    if (strand_type == "unstranded") {
      return "$params.outdir/evidence_data/RNAseq_alignments/STAR/unstranded"
    } else if (strand_type in ["stranded_forward", "stranded_reverse"]) {
      return "$params.outdir/evidence_data/RNAseq_alignments/STAR/stranded"
    }
  }
  cpus 4

  input:
    val(star_database)
    tuple val(sample_ID), val(library_layout), path(reads), val(strand_type)

  output:
    tuple val(sample_ID), path("${sample_ID}_Aligned.sortedByCoord.out.bam"), val(strand_type), emit: samples_aligned

  script:
    def basename_database = task.ext.prefix ?: "${star_database.getName()}"
    """
    DATE=\$(date "+%Y-%m-%d %H:%M:%S")
    echo "[\$DATE] Running STAR alignment on $sample_ID"

    if [[ $library_layout == "paired" ]]
    then

      CMD="STAR --readFilesCommand zcat --genomeDir /star_databases/${basename_database} --runThreadN ${task.cpus} --readFilesIn ${reads[0]} ${reads[1]} --outFileNamePrefix ${sample_ID}_ --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard --outSAMstrandField intronMotif"
      echo "[\$DATE] Executing: \$CMD"

      STAR --readFilesCommand zcat --genomeDir /star_databases/${basename_database} --runThreadN ${task.cpus} --readFilesIn ${reads[0]} ${reads[1]} --outFileNamePrefix ${sample_ID}_ --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard --outSAMstrandField intronMotif
    elif [[ $library_layout == "single" ]]
    then

      CMD="STAR --readFilesCommand zcat --genomeDir /star_databases/${basename_database} --runThreadN ${task.cpus} --readFilesIn ${reads} --outFileNamePrefix ${sample_ID}_ --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard --outSAMstrandField intronMotif"
      echo "[\$DATE] Executing: \$CMD"

      STAR --readFilesCommand zcat --genomeDir /star_databases/${basename_database} --runThreadN ${task.cpus} --readFilesIn ${reads} --outFileNamePrefix ${sample_ID}_ --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard --outSAMstrandField intronMotif
    fi
    """
}
