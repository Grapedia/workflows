// AUGUSTUS can be run directly using the pipeline BRAKER3
process braker3_prediction_with_long_reads {
  label 'process_prediction'

  tag "Executing BRAKER3/AUGUSTUS-Genemark prediction"
  container params.container_braker3
  publishDir "${params.output_dir}", mode: 'copy', saveAs: { filename ->
    if (filename in ['augustus.hints.gff3', 'genemark.gtf', 'genemark_supported.gtf', 'braker.gff3', 'versions.yml']) {
      return filename
    }
    return null
  }
  publishDir "${params.output_dir}/intermediate_files/braker3", mode: 'copy', enabled: params.publish_intermediates, saveAs: { filename ->
    if (filename in ['braker3_run.log', 'braker3_command.txt', 'braker3_inputs.tsv', 'braker.log', 'errors']) {
      return filename
    }
    return null
  }

  input:
    path(genome)
    path(protein_fastas, stageAs: "protein_fastas/*")
    path(bam_short, stageAs: "bam_short/*")
    path(bam_long, stageAs: "bam_long/*")
    path(braker3_runner)

  output:
    path "augustus.hints.gff3", emit: augustus_gff
    path "genemark.gtf", emit: genemark_gtf
    path "genemark_supported.gtf", emit: genemark_supported_gtf
    path "braker.gff3", emit: braker_gff
    path "braker3_run.log", emit: debug_log
    path "braker3_command.txt", emit: command_log
    path "braker3_inputs.tsv", emit: input_manifest
    path "braker.log", optional: true, emit: braker_log
    path "errors", optional: true, emit: errors
    path "versions.yml", emit: versions



  script:
    """
    set -euo pipefail
    bash "${braker3_runner}" "${task.process}" "${genome}" "${task.cpus}" "${task.memory ?: 'not_set'}" "true"
    """

  stub:
    """
    printf "##gff-version 3\\n" > augustus.hints.gff3
    : > genemark.gtf
    cp genemark.gtf genemark_supported.gtf
    printf "##gff-version 3\\n" > braker.gff3
    printf "Stub BRAKER3 run with long reads\\n" > braker3_run.log
    printf "stub\\n" > braker3_command.txt
    printf "type\\tpath\\tsize_bytes\\n" > braker3_inputs.tsv
    printf '"%s":\n  braker: "stub"\n' "${task.process}" > versions.yml
    """
}
