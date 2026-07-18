process agat_convert_gff3_to_gtf {
  label 'process_low'

  tag "Converting ${gff3.name} to GTF for FLAIR"
  container params.container_agat
  publishDir "${params.output_dir}/intermediate_files/liftoff/gff3_to_gtf", mode: "copy", enabled: params.publish_intermediates

  input:
    path(gff3)

  output:
    path "${gff3.simpleName}.gtf", emit: gtf
    path "versions.yml", emit: versions

  script:
    """
    set -euo pipefail
    # A handful of free-text Apollo curator notes have been observed leaking into
    # liftoff_previous_annotations.gff3 as malformed (non-9-column) lines; AGAT's
    # strict GFF3 parser aborts the whole conversion on the first one found, so
    # drop them defensively before use (same guard as modules/mikado.nf).
    awk -F'\\t' 'NF == 9 || /^#/' "${gff3}" > sanitized.gff3
    agat_convert_sp_gff2gtf.pl --gff sanitized.gff3 -o "${gff3.simpleName}.gtf"
    printf '"%s":\n  container: "%s"\n' "${task.process}" "${task.container}" > versions.yml
    """

  stub:
    """
    set -euo pipefail
    printf 'chrStub\\tAGAT\\texon\\t1\\t100\\t.\\t+\\t.\\tgene_id "stub_gene"; transcript_id "stub_tx";\\n' > "${gff3.simpleName}.gtf"
    printf '"%s":\n  container: "%s"\n' "${task.process}" "${task.container}" > versions.yml
    """
}
