process soft_mask_genome {
  label 'process_low'

  tag "Soft-mask the assembly with the EDTA repeat mask"
  container params.container_python
  publishDir "${params.output_dir}/04_evidence", mode: 'copy', saveAs: { filename ->
    filename == 'assembly_softmasked.EDTA.fasta' ? filename : null
  }
  publishDir "${params.output_dir}/05_run_info/intermediate_files/evidence_data/EDTA", mode: 'copy', enabled: params.publish_intermediates, saveAs: { filename ->
    if (filename == 'soft_mask_summary.json') {
      return filename
    }
    return filename == 'versions.yml' ? 'soft_mask_versions.yml' : null
  }

  input:
    path(genome_fasta)
    path(edta_masked_genome)
    path(soft_mask_script)

  output:
    path("assembly_softmasked.EDTA.fasta"), emit: softmasked_genome
    path("soft_mask_summary.json"), emit: summary
    path("versions.yml"), emit: versions

  script:
    """
    set -euo pipefail

    python3 ${soft_mask_script} \\
      --genome ${genome_fasta} \\
      --masked ${edta_masked_genome} \\
      --out assembly_softmasked.EDTA.fasta \\
      --summary soft_mask_summary.json

    python3 --version 2>&1 | sed 's/^/  python: "/; s/\$/"/' | {
      printf '"%s":\\n' "${task.process}"
      cat
    } > versions.yml
    """

  stub:
    """
    set -euo pipefail

    cp ${genome_fasta} assembly_softmasked.EDTA.fasta
    printf '{"stub": true}\\n' > soft_mask_summary.json
    printf '"%s":\\n  python: "stub"\\n' "${task.process}" > versions.yml
    """
}
