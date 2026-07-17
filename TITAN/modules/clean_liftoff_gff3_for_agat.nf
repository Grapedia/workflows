process clean_liftoff_gff3_for_agat {
  label 'process_low'

  tag "Clean Liftoff GFF3 for AGAT"
  container params.container_python
  publishDir "${params.output_dir}/intermediate_files/liftoff/clean_gff3_for_agat", mode: "copy", enabled: params.publish_intermediates

  input:
    path(liftoff_gff3)
    path(clean_liftoff_gff3_script)

  output:
    path("cleaned.OK.gff3"), emit: cleaned_gff3
    path("removed_feature_ids.txt"), emit: removed_feature_ids
    path("versions.yml"), emit: versions

  script:
    """
    set -euo pipefail

    python3 ${clean_liftoff_gff3_script} \\
      --input ${liftoff_gff3} \\
      --output cleaned.OK.gff3 \\
      --removed-ids removed_feature_ids.txt

    python3 --version 2>&1 | sed 's/^/  python: "/; s/\$/"/' | {
      printf '"%s":\\n' "${task.process}"
      cat
    } > versions.yml
    """

  stub:
    """
    set -euo pipefail

    printf "##gff-version 3\\n" > cleaned.OK.gff3
    : > removed_feature_ids.txt
    printf '"%s":\\n  python: "stub"\\n' "${task.process}" > versions.yml
    """
}
