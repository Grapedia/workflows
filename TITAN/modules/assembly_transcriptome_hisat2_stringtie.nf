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
    DEFAULT_GTF="${sample_ID}_transcriptome.gtf"
    ALT_GTF="${sample_ID}_transcriptome.AltCommands.gtf"
    echo "[\$DATE] Running transcriptome assembly with StringTie on HISAT2 $sample_ID"
    bash ${stringtie_transcriptome_script} \\
      ${stringtie_script} \\
      ${stringtie_alt_script} \\
      ${task.cpus} \\
      ${bam_file} \\
      short \\
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
    printf "chr1\\tStringTie\\ttranscript\\t1\\t10\\t.\\t+\\t.\\tgene_id \\"${sample_ID}_hisat_gene\\"; transcript_id \\"${sample_ID}_hisat_tx\\";\\n" > ${sample_ID}_transcriptome.gtf
    printf "chr1\\tStringTie\\ttranscript\\t1\\t10\\t.\\t+\\t.\\tgene_id \\"${sample_ID}_hisat_gene_alt\\"; transcript_id \\"${sample_ID}_hisat_tx_alt\\";\\n" > ${sample_ID}_transcriptome.AltCommands.gtf
    printf '"%s":\n  stringtie: "stub"\n  stringtie_script_sha256: "stub"\n  stringtie_alt_script_sha256: "stub"\n  stringtie_transcriptome_script_sha256: "stub"\n' "${task.process}" > versions.yml
    """
}
