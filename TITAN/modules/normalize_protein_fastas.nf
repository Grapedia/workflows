process normalize_protein_fastas {
  label 'process_low'
  tag "Normalize protein FASTA for ${organism}"
  container params.container_python

  input:
    tuple val(organism), path(protein_fasta)

  output:
    tuple val(organism), path("protein_*.fa"), emit: normalized_fastas
    path "versions.yml", emit: versions

  script:
    """
    set -euo pipefail
    test -s "${protein_fasta}"
    safe_name=\$(printf '%s' "${organism}" | tr -c 'A-Za-z0-9_.-' '_')
    cp "${protein_fasta}" "protein_\${safe_name}.fa"
    printf '"%s":\n  container: "not_recorded"\n' "${task.process}" > versions.yml
    """

  stub:
    """
    set -euo pipefail
    safe_name=\$(printf '%s' "${organism}" | tr -c 'A-Za-z0-9_.-' '_')
    printf ">%s_stub\\nM\\n" "\${safe_name}" > "protein_\${safe_name}.fa"
    printf '"%s":\n  container: "not_recorded"\n' "${task.process}" > versions.yml
    """
}
