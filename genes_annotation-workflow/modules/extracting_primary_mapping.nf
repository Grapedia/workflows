// 3. Extracting primary mapping (not 0x100 flag) from SAM files and retrieve filtered data in BAM format
process extracting_primary_mapping {

    tag "SAM parsing on $sample_ID"
    container 'quay.io/biocontainers/samtools:1.9--h10a08f8_12'
    publishDir "$params.outdir/evidence_data/RNAseq_$stranded_or_unstranded/alignments/new_assembly"
    cpus 4

    input:
      tuple val(sample_ID), val(stranded_or_unstranded), val(paired_or_single), path(sam)

    output:
      tuple val(sample_ID), val(stranded_or_unstranded), val(paired_or_single), file("${sample_ID}_vs_new_assembly.bam")

    script:
    """
    mkdir -p ${workDir}/tmp
    samtools view -F 0x100 -b $sam | samtools sort -o ${sample_ID}_vs_new_assembly.bam -T ${workDir}/tmp/ -
    """
}
