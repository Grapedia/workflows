// Trimming RNAseq data with fastp (much faster than Trimmomatic)
process trimming_fastq {

    tag "FASTP on $sample_ID"
    container 'quay.io/biocontainers/fastp:0.23.2--hb7a2d85_2'
    containerOptions "--volume ${projectDir}/data/RNAseq_data:/RNAseq_data"
    publishDir "$params.outdir/evidence_data/RNAseq_data/trimmed_data"
    cpus 4

    input:
      tuple val(sample_ID), val(SRA_or_FASTQ), val(paired_or_single)

    output:
      tuple val(sample_ID), val(paired_or_single), file("*.trimmed.fastq.gz")

    script:
    """
    if [[ $paired_or_single == "paired" ]]
    then
      fastp --thread ${task.cpus} -i /RNAseq_data/${sample_ID}_1.fastq.gz -I /RNAseq_data/${sample_ID}_2.fastq.gz \
      -o ${sample_ID}_1.trimmed.fastq.gz -O ${sample_ID}_2.trimmed.fastq.gz
    elif [[ $paired_or_single == "single" ]]
    then
      fastp --thread ${task.cpus} -i /RNAseq_data/${sample_ID}.fastq.gz \
      -o ${sample_ID}.trimmed.fastq.gz
    else
      echo "WARNING : \$paired_or_single is not equal to paired or single !"
    fi
    """
}
