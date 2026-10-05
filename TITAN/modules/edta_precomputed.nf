// Stands in for the EDTA process when --edta_precomputed_dir is set: EDTA takes days, so a run
// whose EDTA outputs already exist (e.g. published by an earlier TITAN run) can reuse them.
// The three files are passed through unchanged; only a versions.yml is written.
process edta_precomputed {
  label 'process_low'

  tag "Reuse precomputed EDTA outputs"
  publishDir "${params.output_dir}", mode: 'copy', saveAs: { filename ->
    if (filename == 'assembly_masked.EDTA.fasta') {
      return '04_evidence/assembly_masked.EDTA.fasta'
    }
    if (params.publish_intermediates && filename in ['edta.TElib.fa', 'edta.TEanno.gff3', 'versions.yml']) {
      return "05_run_info/intermediate_files/evidence_data/EDTA/${filename}"
    }
    return null
  }

  input:
    path(masked_genome)
    path(te_annotation)
    path(te_library)

  output:
    path("assembly_masked.EDTA.fasta"), emit: masked_genome
    path("edta.TEanno.gff3"), emit: TE_annotations_gff3
    path("edta.TElib.fa"), emit: TElib_fasta
    path("versions.yml"), emit: versions

  script:
    """
    set -euo pipefail

    printf '"%s":\\n  container: "not_run"\\n  edta_source: "precomputed"\\n  masked_genome: "%s"\\n  te_annotation: "%s"\\n  te_library: "%s"\\n' \\
      "${task.process}" "${masked_genome.toRealPath()}" "${te_annotation.toRealPath()}" "${te_library.toRealPath()}" > versions.yml
    """

  stub:
    """
    set -euo pipefail

    printf '"%s":\\n  container: "not_run"\\n  edta_source: "precomputed"\\n' "${task.process}" > versions.yml
    """
}
