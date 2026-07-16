// 4. Transcriptome assembly with StringTie
process assembly_transcriptome_star_stringtie {
  label 'process_transcriptome'

  tag "STAR/StringTie - short reads on ${sample_ID}"
  container params.container_stringtie
  publishDir { 
    if (strand_type == "unstranded") {
      return "${params.output_dir}/intermediate_files/evidence_data/transcriptomes/StringTie/short_reads/STAR/unstranded"
    } else if (strand_type in ["stranded_forward", "stranded_reverse"]) {
      return "${params.output_dir}/intermediate_files/evidence_data/transcriptomes/StringTie/short_reads/STAR/stranded"
    }
  }
  input:
    tuple val(sample_ID), path(bam_file), val(strand_type)
    path(stringtie_script)
    path(stringtie_alt_script)

  output:
    tuple val(sample_ID), path("${sample_ID}_transcriptome.gtf"), path("${sample_ID}_transcriptome.AltCommands.gtf"), val(strand_type), emit: star_stringtie_transcriptomes

  script:
    """
    DATE=\$(date "+%Y-%m-%d %H:%M:%S")
    echo "[\$DATE] Running STAR/StringTie transcriptome assembly on $sample_ID"
    CMD="${stringtie_script} -t ${task.cpus} -o ${sample_ID}_transcriptome.gtf -b ${bam_file} -r short"
    echo "[\$DATE] Executing: \$CMD"
    ${stringtie_script} -t ${task.cpus} -o ${sample_ID}_transcriptome.gtf -b ${bam_file} -r short
    CMD="${stringtie_alt_script} -t ${task.cpus} -o ${sample_ID}_transcriptome.AltCommands.gtf -b ${bam_file} -r short"
    echo "[\$DATE] Executing: \$CMD"
    ${stringtie_alt_script} -t ${task.cpus} -o ${sample_ID}_transcriptome.AltCommands.gtf -b ${bam_file} -r short
    """

  stub:
    """
    printf "chr1\\tStringTie\\ttranscript\\t1\\t10\\t.\\t+\\t.\\tgene_id \\"${sample_ID}_gene\\"; transcript_id \\"${sample_ID}_tx\\";\\n" > ${sample_ID}_transcriptome.gtf
    printf "chr1\\tStringTie\\ttranscript\\t1\\t10\\t.\\t+\\t.\\tgene_id \\"${sample_ID}_gene_alt\\"; transcript_id \\"${sample_ID}_tx_alt\\";\\n" > ${sample_ID}_transcriptome.AltCommands.gtf
    """
}
