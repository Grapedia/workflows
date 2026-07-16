// 4. Transcriptome assembly with StringTie
process assembly_transcriptome_hisat2_stringtie {
  label 'process_transcriptome'

  tag "Hisat2/StringTie - short reads on ${sample_ID}"
  container params.container_stringtie
  publishDir {
    if (strand_type == "unstranded") {
      return "${params.output_dir}/intermediate_files/evidence_data/transcriptomes/StringTie/short_reads/HISAT2/unstranded"
    } else if (strand_type in ["stranded_forward", "stranded_reverse"]) {
      return "${params.output_dir}/intermediate_files/evidence_data/transcriptomes/StringTie/short_reads/HISAT2/stranded"
    }
  }, mode: "copy", enabled: params.publish_intermediates
  input:
    tuple val(sample_ID), path(bam_file), val(strand_type)
    path(stringtie_script)
    path(stringtie_alt_script)
    path(stringtie_transcriptome_script)

  output:
    tuple val(sample_ID), path("${sample_ID}_transcriptome.gtf"), path("${sample_ID}_transcriptome.AltCommands.gtf"), val(strand_type), emit: hisat2_stringtie_transcriptomes
    path "versions.yml", emit: versions



  script:
    """
    set -euo pipefail
    DATE=\$(date "+%Y-%m-%d %H:%M:%S")
    echo "[\$DATE] Running transcriptome assembly with StringTie on HISAT2 $sample_ID"
    bash ${stringtie_transcriptome_script} \\
      ${stringtie_script} \\
      ${stringtie_alt_script} \\
      ${task.cpus} \\
      ${bam_file} \\
      short \\
      ${sample_ID}_transcriptome.gtf \\
      ${sample_ID}_transcriptome.AltCommands.gtf
    printf '"%s":\n  container: "not_recorded"\n' "${task.process}" > versions.yml
    """

  stub:
    """
    printf "chr1\\tStringTie\\ttranscript\\t1\\t10\\t.\\t+\\t.\\tgene_id \\"${sample_ID}_hisat_gene\\"; transcript_id \\"${sample_ID}_hisat_tx\\";\\n" > ${sample_ID}_transcriptome.gtf
    printf "chr1\\tStringTie\\ttranscript\\t1\\t10\\t.\\t+\\t.\\tgene_id \\"${sample_ID}_hisat_gene_alt\\"; transcript_id \\"${sample_ID}_hisat_tx_alt\\";\\n" > ${sample_ID}_transcriptome.AltCommands.gtf
    printf '"%s":\n  container: "not_recorded"\n' "${task.process}" > versions.yml
    """
}
