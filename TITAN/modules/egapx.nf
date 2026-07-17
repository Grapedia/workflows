process egapx {
  label 'process_prediction'

  tag "Executing NCBI egapx gene annotation pipeline ..."
  // Intentionally no process container: this wrapper launches the official EGAPx
  // nested Nextflow runner, which must access the host Docker/Apptainer runtime.
  publishDir "${params.output_dir}/egapx", mode: 'copy'
  input:
    path egapx_paramfile

  output:
    path "egapx.complete.genomic.gff3", emit: gff3
    path "egapx.complete.genomic.gtf", emit: gtf
    path "egapx.complete.proteins.faa", emit: proteins
    path "egapx.complete.cds.fna", emit: cds
    path "egapx.complete.transcripts.fna", emit: transcripts
    path "egapx.annotated_genome.asn", emit: annotated_genome_asn
    path "egapx_out", emit: output_dir
    path "versions.yml", emit: versions

  script:
    """
    set -euo pipefail
    DATE=\$(date "+%Y-%m-%d %H:%M:%S")

    EGAPX_VERSION="${params.egapx_version}"
    EGAPX_REVISION="${params.egapx_revision}"
    EGAPX_CONTAINER="${params.container_egapx}"
    EGAPX_EXECUTOR="${params.egapx_executor}"
    EGAPX_PYTHON="${params.egapx_python}"
    EGAPX_DATA_VERSION="${params.egapx_data_version}"
    EGAPX_RUNNER_DIR="${params.egapx_runner_dir}"
    EGAPX_LOCAL_CACHE_DIR="${params.egapx_local_cache_dir}"
    EGAPX_CONFIG_DIR="${params.egapx_config_dir}"

    for required_command in "\$EGAPX_PYTHON" curl tar "\$EGAPX_EXECUTOR"; do
      if ! command -v "\$required_command" >/dev/null 2>&1; then
        echo "ERROR: EGAPx wrapper requires '\$required_command' on the host PATH because it launches a nested Nextflow workflow." >&2
        exit 1
      fi
    done

    mkdir -p egapx_runner egapx_work egapx_out

    if [[ -n "\$EGAPX_RUNNER_DIR" && "\$EGAPX_RUNNER_DIR" != "false" ]]; then
      echo "[\$DATE] Using pre-staged EGAPx runner from \$EGAPX_RUNNER_DIR"
      cp -R "\$EGAPX_RUNNER_DIR"/. egapx_runner/
    else
      echo "[\$DATE] Downloading EGAPx runner revision \$EGAPX_REVISION. For strict offline reproducibility, provide --egapx_runner_dir."
      curl -fsSL "https://github.com/ncbi/egapx/archive/refs/tags/\${EGAPX_REVISION}.tar.gz" \\
        | tar -xz --strip-components=1 -C egapx_runner
    fi

    if [[ -n "\$EGAPX_LOCAL_CACHE_DIR" && "\$EGAPX_LOCAL_CACHE_DIR" != "false" ]]; then
      mkdir -p "\$EGAPX_LOCAL_CACHE_DIR"
      local_cache="\$EGAPX_LOCAL_CACHE_DIR"
    else
      mkdir -p egapx_local_cache
      local_cache="\$PWD/egapx_local_cache"
    fi

    printf "process.container = '%s'\\n" "\$EGAPX_CONTAINER" > egapx_runner/ui/assets/config/docker_image.config

    CMD=(
      "\$EGAPX_PYTHON" egapx_runner/ui/egapx.py
      ${egapx_paramfile}
      -e "\$EGAPX_EXECUTOR"
      -w "\$PWD/egapx_work"
      -o "\$PWD/egapx_out"
      -lc "\$local_cache"
      -dv "\$EGAPX_DATA_VERSION"
    )
    if [[ -n "\$EGAPX_CONFIG_DIR" && "\$EGAPX_CONFIG_DIR" != "false" ]]; then
      if [[ ! -d "\$EGAPX_CONFIG_DIR" ]]; then
        echo "ERROR: EGAPx config directory does not exist: \$EGAPX_CONFIG_DIR" >&2
        exit 1
      fi
      CMD+=(-c "\$EGAPX_CONFIG_DIR")
    fi
    echo "[\$DATE] Executing: \${CMD[*]}"
    "\${CMD[@]}"

    copy_required_egapx_output() {
      local source_path="\$1"
      local output_path="\$2"
      local label="\$3"
      if [[ ! -s "\$source_path" ]]; then
        echo "ERROR: EGAPx did not produce expected \$label output: \$source_path" >&2
        echo "       Check egapx_out/nextflow logs, executor '\$EGAPX_EXECUTOR', data version '\$EGAPX_DATA_VERSION', and container '\$EGAPX_CONTAINER'." >&2
        exit 1
      fi
      cp "\$source_path" "\$output_path"
    }

    copy_required_egapx_output egapx_out/complete.genomic.gff egapx.complete.genomic.gff3 "complete genomic GFF"
    copy_required_egapx_output egapx_out/complete.genomic.gtf egapx.complete.genomic.gtf "complete genomic GTF"
    copy_required_egapx_output egapx_out/complete.proteins.faa egapx.complete.proteins.faa "protein FASTA"
    copy_required_egapx_output egapx_out/complete.cds.fna egapx.complete.cds.fna "CDS FASTA"
    copy_required_egapx_output egapx_out/complete.transcripts.fna egapx.complete.transcripts.fna "transcript FASTA"
    copy_required_egapx_output egapx_out/annotated_genome.asn egapx.annotated_genome.asn "annotated genome ASN"

    printf '"%s":\\n  egapx: "%s"\\n  egapx_runner_revision: "%s"\\n  egapx_container: "%s"\\n  egapx_data_version: "%s"\\n' \\
      "${task.process}" "\$EGAPX_VERSION" "\$EGAPX_REVISION" "\$EGAPX_CONTAINER" "\$EGAPX_DATA_VERSION" > versions.yml
    """

  stub:
    """
    mkdir -p egapx_out/GNOMON egapx_out/stats egapx_out/validated egapx_out/nextflow
    printf "##gff-version 3\\nchr1\\tEGAPx\\tgene\\t1\\t10\\t.\\t+\\t.\\tID=egapx_stub_gene\\n" > egapx.complete.genomic.gff3
    printf "chr1\\tEGAPx\\tgene\\t1\\t10\\t.\\t+\\t.\\tgene_id \\"egapx_stub_gene\\";\\n" > egapx.complete.genomic.gtf
    printf ">egapx_stub_protein\\nM\\n" > egapx.complete.proteins.faa
    printf ">egapx_stub_cds\\nATG\\n" > egapx.complete.cds.fna
    printf ">egapx_stub_transcript\\nATG\\n" > egapx.complete.transcripts.fna
    printf "EGAPx ASN stub\\n" > egapx.annotated_genome.asn
    cp egapx.complete.genomic.gff3 egapx_out/complete.genomic.gff
    cp egapx.complete.genomic.gtf egapx_out/complete.genomic.gtf
    cp egapx.complete.proteins.faa egapx_out/complete.proteins.faa
    cp egapx.complete.cds.fna egapx_out/complete.cds.fna
    cp egapx.complete.transcripts.fna egapx_out/complete.transcripts.fna
    cp egapx.annotated_genome.asn egapx_out/annotated_genome.asn
    printf '"%s":\\n  egapx: "%s"\\n  egapx_runner_revision: "%s"\\n  egapx_container: "%s"\\n  egapx_data_version: "%s"\\n' \\
      "${task.process}" "${params.egapx_version}" "${params.egapx_revision}" "${params.container_egapx}" "${params.egapx_data_version}" > versions.yml
    """
}
