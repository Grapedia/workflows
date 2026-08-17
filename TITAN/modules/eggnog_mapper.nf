process eggnog_mapper {
  label 'process_aegis'

  tag "Executing eggnog-mapper on $proteins_file_all and $proteins_file_main"
  container params.container_eggnog_mapper
  publishDir "${params.output_dir}/02_functional_annotation/eggnog", mode: 'copy', saveAs: { filename ->
    if (filename in [
      'final_annotation_proteins_all.emapper.annotations',
      'final_annotation_proteins_main.emapper.annotations',
      'final_annotation_proteins_all.emapper.seed_orthologs',
      'final_annotation_proteins_main.emapper.seed_orthologs',
      'final_annotation_proteins_all.emapper.orthologs',
      'final_annotation_proteins_main.emapper.orthologs',
      'final_annotation_proteins_all.emapper.annotations.xlsx',
      'final_annotation_proteins_main.emapper.annotations.xlsx',
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
    path "final_annotation_proteins_all.emapper.annotations.xlsx", emit: proteins_all_xlsx
    path "final_annotation_proteins_main.emapper.annotations.xlsx", emit: proteins_main_xlsx
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
        touch final_annotation_proteins_all.emapper.annotations.xlsx
        touch final_annotation_proteins_main.emapper.annotations.xlsx
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

    # eggNOG-mapper keys every output (annotations/seed_orthologs/orthologs
    # query column) on the input protein FASTA's own header
    # (<gene_id>_t<NNN>_CDS<N>.prot, a CDS-record ID), not the gene ID
    # downstream analyses join on. Collapse it back to the bare gene ID
    # everywhere it appears (comment/header lines never match, so they pass
    # through untouched), then regenerate the .xlsx exports from the
    # now-cleaned annotations TSVs so they stay in sync.
    sed -i -E 's/(Vitvi[A-Za-z0-9]*g[0-9]+)_t[0-9]+_CDS[0-9]+\\.prot/\\1/g' \\
      final_annotation_proteins_all.emapper.annotations \\
      final_annotation_proteins_main.emapper.annotations \\
      final_annotation_proteins_all.emapper.seed_orthologs \\
      final_annotation_proteins_main.emapper.seed_orthologs \\
      final_annotation_proteins_all.emapper.orthologs \\
      final_annotation_proteins_main.emapper.orthologs

    python3 - <<'PY'
def regenerate_xlsx(tsv_path, xlsx_path):
    rows = [line.rstrip("\n").split("\t") for line in open(tsv_path, encoding="utf-8")]
    try:
        import xlsxwriter
        workbook = xlsxwriter.Workbook(xlsx_path)
        sheet = workbook.add_worksheet()
        for r, row in enumerate(rows):
            for c, value in enumerate(row):
                sheet.write(r, c, value)
        workbook.close()
    except ImportError:
        import openpyxl
        workbook = openpyxl.Workbook()
        sheet = workbook.active
        for row in rows:
            sheet.append(row)
        workbook.save(xlsx_path)

regenerate_xlsx("final_annotation_proteins_all.emapper.annotations", "final_annotation_proteins_all.emapper.annotations.xlsx")
regenerate_xlsx("final_annotation_proteins_main.emapper.annotations", "final_annotation_proteins_main.emapper.annotations.xlsx")
PY

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
    touch final_annotation_proteins_all.emapper.annotations.xlsx
    touch final_annotation_proteins_main.emapper.annotations.xlsx
    printf '"%s":\n  eggnog_mapper: "stub"\n  container: "%s"\n' "${task.process}" "${task.container}" > versions.yml
    """
}
