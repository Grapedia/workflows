// 4. Transcriptome assembly with StringTie
process assembly_transcriptome_hisat2_stringtie {

  tag "Hisat2/StringTie - short reads on ${sample_ID}"
  container params.container_stringtie
  publishDir { 
    if (strand_type == "unstranded") {
      return "${params.output_dir}/intermediate_files/evidence_data/transcriptomes/StringTie/short_reads/HISAT2/unstranded"
    } else if (strand_type in ["stranded_forward", "stranded_reverse"]) {
      return "${params.output_dir}/intermediate_files/evidence_data/transcriptomes/StringTie/short_reads/HISAT2/stranded"
    }
  }
  cpus 4

  input:
    tuple val(sample_ID), path(bam_file), val(strand_type)
    path(stringtie_script)
    path(stringtie_alt_script)

  output:
    tuple val(sample_ID), path("${sample_ID}_transcriptome.gtf"), path("${sample_ID}_transcriptome.AltCommands.gtf"), val(strand_type), emit: hisat2_stringtie_transcriptomes

  script:
    """
    DATE=\$(date "+%Y-%m-%d %H:%M:%S")
    echo "[\$DATE] Running transcriptome assembly with StringTie on HISAT2 $sample_ID"
    CMD="${stringtie_script} -t ${task.cpus} -o ${sample_ID}_transcriptome.gtf -b ${bam_file} -r short"
    echo "[\$DATE] Executing: \$CMD"
    ${stringtie_script} -t ${task.cpus} -o ${sample_ID}_transcriptome.gtf -b ${bam_file} -r short
    CMD="${stringtie_alt_script} -t ${task.cpus} -o ${sample_ID}_transcriptome.AltCommands.gtf -b ${bam_file} -r short"
    echo "[\$DATE] Executing: \$CMD"
    ${stringtie_alt_script} -t ${task.cpus} -o ${sample_ID}_transcriptome.AltCommands.gtf -b ${bam_file} -r short
    """

  stub:
    """
    printf "chr1\\tStringTie\\ttranscript\\t1\\t10\\t.\\t+\\t.\\tgene_id \\"${sample_ID}_hisat_gene\\"; transcript_id \\"${sample_ID}_hisat_tx\\";\\n" > ${sample_ID}_transcriptome.gtf
    printf "chr1\\tStringTie\\ttranscript\\t1\\t10\\t.\\t+\\t.\\tgene_id \\"${sample_ID}_hisat_gene_alt\\"; transcript_id \\"${sample_ID}_hisat_tx_alt\\";\\n" > ${sample_ID}_transcriptome.AltCommands.gtf
    """
}
