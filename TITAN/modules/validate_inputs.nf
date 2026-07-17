process validate_inputs {
  label 'process_low'
  tag "TITAN input validation"
  container params.container_python

  input:
    path(new_assembly)
    path(previous_assembly)
    path(previous_annotations)
    path(rnaseq_samplesheet)
    val(rnaseq_data_dir)
    path(protein_samplesheet)
    path(egapx_paramfile)
    val(egapx_paramfile_source)
    val(egapx_executor)
    val(psiclass_vd)
    val(psiclass_c)
    val(egapx_version)
    val(egapx_revision)
    val(egapx_container)
    val(egapx_data_version)
    val(aegis_version)
    val(aegis_container)
    val(run_eggnog_mapper)
    val(eggnog_data_dir)
    val(run_helixer)
    val(helixer_model_dir)
    val(helixer_model)
    val(run_interproscan)
    val(interproscan_data_dir)

  output:
    path "validated_inputs.ok", emit: ok
    path "versions.yml", emit: versions

  script:
    """
    set -euo pipefail

    egapx_version="${egapx_version}"
    egapx_revision="${egapx_revision}"
    egapx_container="${egapx_container}"
    egapx_data_version="${egapx_data_version}"
    aegis_version="${aegis_version}"
    aegis_container="${aegis_container}"

    for name in egapx_version egapx_revision egapx_container egapx_data_version aegis_version aegis_container; do
      value="\${!name}"
      if [[ -z "\${value}" || "\${value}" == "false" || "\${value}" == "true" ]]; then
        echo "ERROR: runtime parameter \${name} must be explicitly configured" >&2
        exit 1
      fi
    done

    python3 ${projectDir}/scripts/validate_inputs.py \\
      --project-dir "${projectDir}" \\
      --new-assembly "${new_assembly}" \\
      --previous-assembly "${previous_assembly}" \\
      --previous-annotations "${previous_annotations}" \\
      --rnaseq-samplesheet "${rnaseq_samplesheet}" \\
      --rnaseq-data-dir "${rnaseq_data_dir}" \\
      --protein-samplesheet "${protein_samplesheet}" \\
      --egapx-paramfile "${egapx_paramfile_source}" \\
      --egapx-executor "${egapx_executor}" \\
      --psiclass-vd "${psiclass_vd}" \\
      --psiclass-c "${psiclass_c}" \\
      --run-eggnog-mapper "${run_eggnog_mapper}" \\
      --eggnog-data-dir "${eggnog_data_dir}" \\
      --run-helixer "${run_helixer}" \\
      --helixer-model-dir "${helixer_model_dir}" \\
      --helixer-model "${helixer_model}" \\
      --run-interproscan "${run_interproscan}" \\
      --interproscan-data-dir "${interproscan_data_dir}"

    touch validated_inputs.ok
    printf '"%s":\n  container: "not_recorded"\n' "${task.process}" > versions.yml
    """

  stub:
    """
    set -euo pipefail

    egapx_version="${egapx_version}"
    egapx_revision="${egapx_revision}"
    egapx_container="${egapx_container}"
    egapx_data_version="${egapx_data_version}"
    aegis_version="${aegis_version}"
    aegis_container="${aegis_container}"

    for name in egapx_version egapx_revision egapx_container egapx_data_version aegis_version aegis_container; do
      value="\${!name}"
      if [[ -z "\${value}" || "\${value}" == "false" || "\${value}" == "true" ]]; then
        echo "ERROR: runtime parameter \${name} must be explicitly configured" >&2
        exit 1
      fi
    done

    python3 ${projectDir}/scripts/validate_inputs.py \\
      --project-dir "${projectDir}" \\
      --new-assembly "${new_assembly}" \\
      --previous-assembly "${previous_assembly}" \\
      --previous-annotations "${previous_annotations}" \\
      --rnaseq-samplesheet "${rnaseq_samplesheet}" \\
      --rnaseq-data-dir "${rnaseq_data_dir}" \\
      --protein-samplesheet "${protein_samplesheet}" \\
      --egapx-paramfile "${egapx_paramfile_source}" \\
      --egapx-executor "${egapx_executor}" \\
      --psiclass-vd "${psiclass_vd}" \\
      --psiclass-c "${psiclass_c}" \\
      --run-eggnog-mapper "${run_eggnog_mapper}" \\
      --eggnog-data-dir "${eggnog_data_dir}" \\
      --run-helixer "${run_helixer}" \\
      --helixer-model-dir "${helixer_model_dir}" \\
      --helixer-model "${helixer_model}" \\
      --run-interproscan "${run_interproscan}" \\
      --interproscan-data-dir "${interproscan_data_dir}"

    touch validated_inputs.ok
    printf '"%s":\n  container: "not_recorded"\n' "${task.process}" > versions.yml
    """
}
