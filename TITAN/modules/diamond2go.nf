process diamond2go {
  label 'process_aegis'

  tag "Executing diamond2go annotation on $proteins_file_all and $proteins_file_main"
  container params.container_diamond2go
  publishDir "${params.output_dir}/Diamond2GO_outputs", mode: 'copy', saveAs: { filename ->
    if (filename in ['final_annotation_proteins_all.diamond2go.tsv', 'final_annotation_proteins_main.diamond2go.tsv', 'versions.yml']) {
      return filename
    }
    return null
  }

  input:
    path(proteins_file_all)
    path(proteins_file_main)

  output:
    path "final_annotation_proteins_all.diamond2go.tsv", emit: proteins_all_diamond
    path "final_annotation_proteins_main.diamond2go.tsv", emit: proteins_main_diamond
    path "versions.yml", emit: versions



  script:
    """
    set -euo pipefail

    D2GO_SCRIPT="/Diamond2GO/Diamond2go.pl"
    D2GO_DB="/Diamond2GO/resources/nr_clean_d2go.dmnd"

    require_file() {
        local label="\$1"
        local file="\$2"
        if [[ ! -s "\$file" ]]; then
            echo "ERROR: missing or empty \${label}: \${file}" >&2
            exit 1
        fi
    }

    select_diamond2go_output() {
        local run_dir="\$1"
        local label="\$2"
        local selected=""

        selected=\$(find "\$run_dir" -maxdepth 1 -type f \\( -name '*diamond*processed*' -o -name '*diamond2go*' -o -name '*-diamond*' -o -name '*.diamond.*' \\) | sort | head -n 1)
        if [[ -z "\$selected" || ! -s "\$selected" ]]; then
            echo "ERROR: Diamond2GO did not produce a non-empty \${label} output in \${run_dir}" >&2
            find "\$run_dir" -maxdepth 1 -type f | sort >&2
            exit 1
        fi
        printf '%s\\n' "\$selected"
    }

    run_diamond2go() {
        local label="\$1"
        local input_fasta="\$2"
        local output_tsv="\$3"
        local run_dir="diamond2go_\${label}"

        mkdir -p "\$run_dir"
        cp "\$input_fasta" "\$run_dir/query.fasta"
        pushd "\$run_dir" >/dev/null
        perl "\$D2GO_SCRIPT" -d "\$D2GO_DB" -q query.fasta -t protein -r "${task.cpus}"
        popd >/dev/null

        selected=\$(select_diamond2go_output "\$run_dir" "\$label")
        cp "\$selected" "\$output_tsv"
        test -s "\$output_tsv"
    }

    require_file "Diamond2GO script" "\$D2GO_SCRIPT"
    require_file "Diamond2GO database" "\$D2GO_DB"
    require_file "all proteins FASTA" "${proteins_file_all}"
    require_file "main proteins FASTA" "${proteins_file_main}"

    run_diamond2go all "${proteins_file_all}" final_annotation_proteins_all.diamond2go.tsv
    run_diamond2go main "${proteins_file_main}" final_annotation_proteins_main.diamond2go.tsv

    diamond_version=\$(diamond version 2>/dev/null | head -n 1 || true)
    printf '"%s":\n  diamond2go: "%s"\n  diamond: "%s"\n  threads: "%s"\n  container: "%s"\n' \\
      "${task.process}" "Diamond2go.pl" "\${diamond_version:-unknown}" "${task.cpus}" "${task.container}" > versions.yml
    """

  stub:
    """
    set -euo pipefail
    printf "query\\tsubject\\n" > final_annotation_proteins_all.diamond2go.tsv
    printf "query\\tsubject\\n" > final_annotation_proteins_main.diamond2go.tsv
    printf '"%s":\n  diamond2go: "stub"\n  threads: "%s"\n' "${task.process}" "${task.cpus}" > versions.yml
    """
}
