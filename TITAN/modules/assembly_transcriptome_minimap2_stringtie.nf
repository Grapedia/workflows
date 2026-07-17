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
    DEFAULT_GTF="${sample_ID}_transcriptome.gtf"
    ALT_GTF="${sample_ID}_transcriptome.AltCommands.gtf"
    echo "[\$DATE] Running Minimap2/StringTie - long reads transcriptome assembly on $sample_ID"
    bash ${stringtie_transcriptome_script} \\
      ${stringtie_script} \\
      ${stringtie_alt_script} \\
      ${task.cpus} \\
      ${bam_file} \\
      long \\
      "\${DEFAULT_GTF}" \\
      "\${ALT_GTF}"
    {
      printf '"%s":\n' "${task.process}"
      stringtie --version 2>&1 | awk '{ printf "  stringtie: \\"%s\\"\\n", \$0 }'
      sha256sum ${stringtie_script} | awk '{ printf "  stringtie_script_sha256: \\"%s\\"\\n", \$1 }'
      sha256sum ${stringtie_alt_script} | awk '{ printf "  stringtie_alt_script_sha256: \\"%s\\"\\n", \$1 }'
      sha256sum ${stringtie_transcriptome_script} | awk '{ printf "  stringtie_transcriptome_script_sha256: \\"%s\\"\\n", \$1 }'
    } > versions.yml
    """

  stub:
    """
    set -euo pipefail
    printf "chr1\\tStringTie\\ttranscript\\t1\\t10\\t.\\t+\\t.\\tgene_id \\"${sample_ID}_long_gene\\"; transcript_id \\"${sample_ID}_long_tx\\";\\n" > ${sample_ID}_transcriptome.gtf
    printf "chr1\\tStringTie\\ttranscript\\t1\\t10\\t.\\t+\\t.\\tgene_id \\"${sample_ID}_long_gene_alt\\"; transcript_id \\"${sample_ID}_long_tx_alt\\";\\n" > ${sample_ID}_transcriptome.AltCommands.gtf
    printf '"%s":\n  stringtie: "stub"\n  stringtie_script_sha256: "stub"\n  stringtie_alt_script_sha256: "stub"\n  stringtie_transcriptome_script_sha256: "stub"\n' "${task.process}" > versions.yml
    """
}
