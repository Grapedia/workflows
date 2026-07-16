// 4. Transcriptome assembly with StringTie
process assembly_transcriptome_minimap2_stringtie {
  label 'process_transcriptome'

  tag "Minimap2/StringTie - long reads"
  container params.container_stringtie
  publishDir "${params.output_dir}/intermediate_files/evidence_data/transcriptomes/StringTie/long_reads", mode: "copy", enabled: params.publish_intermediates
  input:
    tuple val(sample_ID), path(bam_file)
    path(stringtie_script)
    path(stringtie_alt_script)
    path(stringtie_transcriptome_script)

  output:
    tuple val(sample_ID), path("${sample_ID}_transcriptome.gtf"), path("${sample_ID}_transcriptome.AltCommands.gtf"), emit: minimap2_stringtie_transcriptomes
    path "versions.yml", emit: versions



  script:
    """
    set -euo pipefail
    DATE=\$(date "+%Y-%m-%d %H:%M:%S")
    echo "[\$DATE] Running Minimap2/StringTie - long reads transcriptome assembly on $sample_ID"
    bash ${stringtie_transcriptome_script} \\
      ${stringtie_script} \\
      ${stringtie_alt_script} \\
      ${task.cpus} \\
      ${bam_file} \\
      long \\
      ${sample_ID}_transcriptome.gtf \\
      ${sample_ID}_transcriptome.AltCommands.gtf
    printf '"%s":\n  container: "not_recorded"\n' "${task.process}" > versions.yml
    """

  stub:
    """
    printf "chr1\\tStringTie\\ttranscript\\t1\\t10\\t.\\t+\\t.\\tgene_id \\"${sample_ID}_long_gene\\"; transcript_id \\"${sample_ID}_long_tx\\";\\n" > ${sample_ID}_transcriptome.gtf
    printf "chr1\\tStringTie\\ttranscript\\t1\\t10\\t.\\t+\\t.\\tgene_id \\"${sample_ID}_long_gene_alt\\"; transcript_id \\"${sample_ID}_long_tx_alt\\";\\n" > ${sample_ID}_transcriptome.AltCommands.gtf
    printf '"%s":\n  container: "not_recorded"\n' "${task.process}" > versions.yml
    """
}
