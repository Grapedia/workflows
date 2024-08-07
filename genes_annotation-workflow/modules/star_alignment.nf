// 2. Aligning RNAseq data on reference genome with STAR
process star_alignment {

  tag "STAR on ${sample_ID}"
  container 'quay.io/biocontainers/star:2.7.11b--h43eeafb_2'
  containerOptions "--volume $params.outdir/evidence_data/star_databases/:/star_databases"
  publishDir "$params.outdir/evidence_data/RNAseq_$stranded_or_unstranded/alignments/STAR_new_assembly"
  cpus 4

  input:
    val(star_database)
    tuple val(sample_ID), val(stranded_or_unstranded), val(paired_or_single), path(reads)

  output:
    tuple val(sample_ID), val(stranded_or_unstranded), val(paired_or_single), file("${sample_ID}_vs_new_assembly.sam")

  script:
    def basename_database = task.ext.prefix ?: "${star_database.getName()}"
    """
    if [[ $paired_or_single == "paired" ]]
    then
      STAR --readFilesCommand zcat --genomeDir /star_databases/${basename_database} --runThreadN ${task.cpus} --readFilesIn ${reads[0]} ${reads[1]} --outFileNamePrefix ${sample_ID}_ --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard
    elif [[ $paired_or_single == "single" ]]
    then
      STAR --readFilesCommand zcat --genomeDir /star_databases/${basename_database} --runThreadN ${task.cpus} --readFilesIn ${reads} --outFileNamePrefix ${sample_ID}_ --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard
    fi
    """
}
