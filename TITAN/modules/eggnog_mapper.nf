process eggnog_mapper {
  label 'process_aegis'

  tag "Executing eggnog-mapper on $proteins_file_all and $proteins_file_main"
  container params.container_eggnog_mapper
  publishDir "${params.output_dir}/EggNOG_outputs", mode: 'copy', saveAs: { filename ->
    if (filename in [
      'final_annotation_proteins_all.emapper.annotations',
      'final_annotation_proteins_main.emapper.annotations',
      'final_annotation_proteins_all.emapper.seed_orthologs',
      'final_annotation_proteins_main.emapper.seed_orthologs',
      'final_annotation_proteins_all.emapper.orthologs',
      'final_annotation_proteins_main.emapper.orthologs',
      'final_annotation_proteins_all.xlsx',
      'final_annotation_proteins_main.xlsx',
      'versions.yml'
    ]) {
      return filename
    }
    return null
  }

  input:
    path(proteins_file_all)
    path(proteins_file_main)

  output:
    path "final_annotation_proteins_all.emapper.annotations", emit: proteins_all_annotations
    path "final_annotation_proteins_main.emapper.annotations", emit: proteins_main_annotations
    path "final_annotation_proteins_all.emapper.seed_orthologs", emit: proteins_all_seed_orthologs
    path "final_annotation_proteins_main.emapper.seed_orthologs", emit: proteins_main_seed_orthologs
    path "final_annotation_proteins_all.emapper.orthologs", emit: proteins_all_orthologs
    path "final_annotation_proteins_main.emapper.orthologs", emit: proteins_main_orthologs
    path "final_annotation_proteins_all.xlsx", emit: proteins_all_xlsx
    path "final_annotation_proteins_main.xlsx", emit: proteins_main_xlsx
    path "versions.yml", emit: versions

  script:
    """
    set -euo pipefail

    require_file() {
        local label="\$1"
        local file="\$2"
        if [[ ! -s "\$file" ]]; then
            echo "ERROR: missing or empty \${label}: \${file}" >&2
            exit 1
        fi
    }

    if [[ "${params.run_eggnog_mapper}" != "true" ]]; then
        touch final_annotation_proteins_all.emapper.annotations
        touch final_annotation_proteins_main.emapper.annotations
        touch final_annotation_proteins_all.emapper.seed_orthologs
        touch final_annotation_proteins_main.emapper.seed_orthologs
        touch final_annotation_proteins_all.emapper.orthologs
        touch final_annotation_proteins_main.emapper.orthologs
        touch final_annotation_proteins_all.xlsx
        touch final_annotation_proteins_main.xlsx
        printf '"%s":\n  eggnog_mapper: "skipped"\n  sensmode: "%s"\n  container: "%s"\n' \
            "${task.process}" "${params.eggnog_mapper_sensmode}" "${task.container}" > versions.yml
        exit 0
    fi

    require_file "all proteins FASTA" "${proteins_file_all}"
    require_file "main proteins FASTA" "${proteins_file_main}"

    if [[ -z "${params.eggnog_data_dir}" || "${params.eggnog_data_dir}" == "false" ]]; then
        echo "ERROR: eggnog_data_dir must be provided when run_eggnog_mapper=true" >&2
        exit 1
    fi

    if [[ -n "${params.eggnog_mapper_tax_scope}" && "${params.eggnog_mapper_tax_scope}" != "false" ]]; then
        emapper.py -i "${proteins_file_all}" --itype proteins -m diamond --cpu ${task.cpus} --data_dir "${params.eggnog_data_dir}" --output final_annotation_proteins_all --output_dir . --excel --report_orthologs --sensmode "${params.eggnog_mapper_sensmode}" --override --tax_scope "${params.eggnog_mapper_tax_scope}"
        emapper.py -i "${proteins_file_main}" --itype proteins -m diamond --cpu ${task.cpus} --data_dir "${params.eggnog_data_dir}" --output final_annotation_proteins_main --output_dir . --excel --report_orthologs --sensmode "${params.eggnog_mapper_sensmode}" --override --tax_scope "${params.eggnog_mapper_tax_scope}"
    else
        emapper.py -i "${proteins_file_all}" --itype proteins -m diamond --cpu ${task.cpus} --data_dir "${params.eggnog_data_dir}" --output final_annotation_proteins_all --output_dir . --excel --report_orthologs --sensmode "${params.eggnog_mapper_sensmode}" --override
        emapper.py -i "${proteins_file_main}" --itype proteins -m diamond --cpu ${task.cpus} --data_dir "${params.eggnog_data_dir}" --output final_annotation_proteins_main --output_dir . --excel --report_orthologs --sensmode "${params.eggnog_mapper_sensmode}" --override
    fi

    test -s final_annotation_proteins_all.emapper.annotations
    test -s final_annotation_proteins_main.emapper.annotations

    printf '"%s":\n  eggnog_mapper: "emapper.py"\n  sensmode: "%s"\n  container: "%s"\n' \
        "${task.process}" "${params.eggnog_mapper_sensmode}" "${task.container}" > versions.yml
    """

  stub:
    """
    set -euo pipefail
    touch final_annotation_proteins_all.emapper.annotations
    touch final_annotation_proteins_main.emapper.annotations
    touch final_annotation_proteins_all.emapper.seed_orthologs
    touch final_annotation_proteins_main.emapper.seed_orthologs
    touch final_annotation_proteins_all.emapper.orthologs
    touch final_annotation_proteins_main.emapper.orthologs
    touch final_annotation_proteins_all.xlsx
    touch final_annotation_proteins_main.xlsx
    printf '"%s":\n  eggnog_mapper: "stub"\n  container: "%s"\n' "${task.process}" "${task.container}" > versions.yml
    """
}
