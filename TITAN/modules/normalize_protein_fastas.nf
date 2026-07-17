process normalize_protein_fastas {
  label 'process_low'
  tag "Normalize protein FASTA for ${organism}"
  container params.container_python

  input:
    tuple val(organism), path(protein_fasta)
    path(clean_protein_script)

  output:
    tuple val(organism), path("protein_*.fa"), emit: normalized_fastas
    path "versions.yml", emit: versions

  script:
    """
    set -euo pipefail
    test -s "${protein_fasta}"
    safe_name=\$(printf '%s' "${organism}" | tr -c 'A-Za-z0-9_.-' '_')
    python3 "${clean_protein_script}" "${protein_fasta}" "protein_\${safe_name}.fa" "\${safe_name}"
    test -s "protein_\${safe_name}.fa"
    printf '"%s":\n  container: "%s"\n  cleaner: "clean_protein_fasta_for_BRAKER3.py"\n' "${task.process}" "${task.container}" > versions.yml
    """

  stub:
    """
    set -euo pipefail
    safe_name=\$(printf '%s' "${organism}" | tr -c 'A-Za-z0-9_.-' '_')
    printf ">%s_stub\\nM\\n" "\${safe_name}" > "protein_\${safe_name}.fa"
    printf '"%s":\n  container: "not_recorded"\n' "${task.process}" > versions.yml
    """
}
