// 4. Transcriptome assembly with StringTie
process assembly_transcriptome_minimap2_stringtie {
  label 'process_transcriptome'

  tag "Minimap2/StringTie - long reads"
  container params.container_stringtie
  publishDir "${params.output_dir}/intermediate_files/evidence_data/transcriptomes/StringTie/long_reads"
  input:
    tuple val(sample_ID), path(bam_file)
    path(stringtie_script)
    path(stringtie_alt_script)

  output:
    tuple val(sample_ID), path("${sample_ID}_transcriptome.gtf"), path("${sample_ID}_transcriptome.AltCommands.gtf"), emit: minimap2_stringtie_transcriptomes
    path "versions.yml", emit: versions



  script:
    """
    set -euo pipefail
    DATE=\$(date "+%Y-%m-%d %H:%M:%S")
    echo "[\$DATE] Running Minimap2/StringTie - long reads transcriptome assembly on $sample_ID"
    CMD="${stringtie_script} -t ${task.cpus} -o ${sample_ID}_transcriptome.gtf -b ${bam_file} -r long"
    echo "[\$DATE] Executing: \$CMD"
    ${stringtie_script} -t ${task.cpus} -o ${sample_ID}_transcriptome.gtf -b ${bam_file} -r long
    CMD="${stringtie_alt_script} -t ${task.cpus} -o ${sample_ID}_transcriptome.AltCommands.gtf -b ${bam_file} -r long"
    echo "[\$DATE] Executing: \$CMD"
    ${stringtie_alt_script} -t ${task.cpus} -o ${sample_ID}_transcriptome.AltCommands.gtf -b ${bam_file} -r long
    printf '"%s":\n  container: "not_recorded"\n' "${task.process}" > versions.yml
    """

  stub:
    """
    printf "chr1\\tStringTie\\ttranscript\\t1\\t10\\t.\\t+\\t.\\tgene_id \\"${sample_ID}_long_gene\\"; transcript_id \\"${sample_ID}_long_tx\\";\\n" > ${sample_ID}_transcriptome.gtf
    printf "chr1\\tStringTie\\ttranscript\\t1\\t10\\t.\\t+\\t.\\tgene_id \\"${sample_ID}_long_gene_alt\\"; transcript_id \\"${sample_ID}_long_tx_alt\\";\\n" > ${sample_ID}_transcriptome.AltCommands.gtf
    printf '"%s":\n  container: "not_recorded"\n' "${task.process}" > versions.yml
    """
}
