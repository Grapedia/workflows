process helixer_prediction {
  label 'process_prediction'

  tag "Executing Helixer ab initio prediction on $masked_genome"
  container params.container_helixer
  publishDir "${params.output_dir}/additional_annotations/helixer", mode: 'copy', saveAs: { filename ->
    if (filename in ['helixer.gff3', 'versions.yml']) {
      return filename
    }
    return null
  }

  input:
    path(masked_genome)

  output:
    path "helixer.gff3", emit: gff3
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

    if [[ "${params.run_helixer}" != "true" ]]; then
        touch helixer.gff3
        printf '"%s":\n  helixer: "skipped"\n  lineage: "%s"\n  container: "%s"\n' \
            "${task.process}" "${params.helixer_model}" "${task.container}" > versions.yml
        exit 0
    fi

    require_file "masked genome" "${masked_genome}"

    if [[ -z "${params.helixer_model_dir}" || "${params.helixer_model_dir}" == "false" ]]; then
        echo "ERROR: helixer_model_dir must be provided when run_helixer=true" >&2
        exit 1
    fi

    if [[ ! -d "${params.helixer_model_dir}/${params.helixer_model}" ]]; then
        echo "ERROR: no pre-fetched Helixer model for lineage '${params.helixer_model}' in ${params.helixer_model_dir}" >&2
        echo "       Run: scripts/download_helixer_model.sh --model-dir ${params.helixer_model_dir} --lineage ${params.helixer_model}" >&2
        exit 1
    fi

    mkdir -p helixer_tmp

    export OMP_NUM_THREADS=${task.cpus}
    export TF_NUM_INTRAOP_THREADS=${task.cpus}
    export TF_NUM_INTEROP_THREADS=1

    Helixer.py \
      --fasta-path "${masked_genome}" \
      --lineage "${params.helixer_model}" \
      --gff-output-path helixer.gff3 \
      --temporary-dir helixer_tmp \
      --downloaded-model-path "${params.helixer_model_dir}"

    test -s helixer.gff3

    printf '"%s":\n  helixer: "Helixer.py"\n  lineage: "%s"\n  gpu: "%s"\n  container: "%s"\n' \
        "${task.process}" "${params.helixer_model}" "${params.helixer_use_gpu}" "${task.container}" > versions.yml
    """

  stub:
    """
    set -euo pipefail
    printf "##gff-version 3\\nchr1\\tHelixer\\tgene\\t1\\t10\\t.\\t+\\t.\\tID=stub_gene\\n" > helixer.gff3
    printf '"%s":\n  helixer: "stub"\n  container: "%s"\n' "${task.process}" "${task.container}" > versions.yml
    """
}
