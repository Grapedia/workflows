process transdecoder_longorfs {
  label 'process_medium'
  tag "TransDecoder LongOrfs on Mikado transcripts"
  container params.container_transdecoder

  input:
    path(prepared_fasta)

  output:
    path "mikado_prepared.fasta.transdecoder_dir", emit: longorfs_dir
    path "mikado_prepared.fasta.transdecoder_dir.__checkpoints_longorfs", optional: true, emit: checkpoints
    path "versions.yml", emit: versions

  script:
    """
    set -euo pipefail
    if [[ "${params.run_mikado}" != "true" || "${params.run_transdecoder}" != "true" || ! -s "${prepared_fasta}" ]]; then
      mkdir -p mikado_prepared.fasta.transdecoder_dir
      printf '"%s":\n  transdecoder_longorfs: "skipped"\n  container: "%s"\n' "${task.process}" "${task.container}" > versions.yml
      exit 0
    fi

    # Use the real scripts; the image's /usr/local/bin symlinks are broken.
    export PATH="/usr/local/opt/transdecoder/util:/usr/local/opt/transdecoder:\${PATH}"

    cp "${prepared_fasta}" mikado_prepared.fasta
    TransDecoder.LongOrfs -t mikado_prepared.fasta

    # TransDecoder.LongOrfs has no --version flag (exits non-zero with
    # "don't understand options" if asked); read the installed conda
    # package version instead.
    transdecoder_version=\$(basename "\$(ls /usr/local/conda-meta/transdecoder-*.json 2>/dev/null | head -n1)" .json | sed -E 's/^transdecoder-([^-]+)-.*/\\1/')
    printf '"%s":\\n  transdecoder_longorfs: "%s"\\n  container: "%s"\\n' \\
      "${task.process}" "\${transdecoder_version:-unknown}" "${task.container}" > versions.yml
    """

  stub:
    """
    set -euo pipefail
    mkdir -p mikado_prepared.fasta.transdecoder_dir
    printf '"%s":\n  transdecoder_longorfs: "stub"\n' "${task.process}" > versions.yml
    """
}

process transdecoder_predict {
  label 'process_medium'
  tag "TransDecoder Predict on Mikado transcripts"
  container params.container_transdecoder
  publishDir "${params.output_dir}/evidence/mikado/transdecoder", mode: 'copy', saveAs: { filename ->
    if (filename in ['mikado_prepared.fasta.transdecoder.bed', 'mikado_prepared.fasta.transdecoder.pep', 'mikado_prepared.fasta.transdecoder.gff3', 'versions.yml']) {
      return filename
    }
    return null
  }

  input:
    path(prepared_fasta)
    path(longorfs_dir, stageAs: "transdecoder_longorfs_input_dir")

  output:
    path "mikado_prepared.fasta.transdecoder.bed", emit: bed
    path "mikado_prepared.fasta.transdecoder.pep", emit: pep
    path "mikado_prepared.fasta.transdecoder.gff3", emit: gff3
    path "versions.yml", emit: versions

  script:
    """
    set -euo pipefail
    if [[ "${params.run_mikado}" != "true" || "${params.run_transdecoder}" != "true" || ! -s "${prepared_fasta}" ]]; then
      : > mikado_prepared.fasta.transdecoder.bed
      printf ">transdecoder_empty\\nM\\n" > mikado_prepared.fasta.transdecoder.pep
      printf "##gff-version 3\\n" > mikado_prepared.fasta.transdecoder.gff3
      printf '"%s":\n  transdecoder_predict: "skipped"\n  container: "%s"\n' "${task.process}" "${task.container}" > versions.yml
      exit 0
    fi

    # Use the real scripts; the image's /usr/local/bin symlinks are broken.
    export PATH="/usr/local/opt/transdecoder/util:/usr/local/opt/transdecoder:\${PATH}"

    cp "${prepared_fasta}" mikado_prepared.fasta
    # -L: longorfs_dir is staged under a different name (see the process
    # input) specifically so this is a real, independent copy. Plain `cp -r`
    # on a symlinked directory recreates a symlink at the destination instead
    # of copying its contents, so TransDecoder.Predict would write its scratch
    # files straight through into transdecoder_longorfs's cached output dir,
    # corrupting it for every future resume.
    cp -rL "${longorfs_dir}" mikado_prepared.fasta.transdecoder_dir
    TransDecoder.Predict -t mikado_prepared.fasta --single_best_only
    test -s mikado_prepared.fasta.transdecoder.bed
    test -s mikado_prepared.fasta.transdecoder.pep
    test -s mikado_prepared.fasta.transdecoder.gff3

    # TransDecoder.Predict has no --version flag either; same conda-meta
    # lookup as transdecoder_longorfs.
    transdecoder_version=\$(basename "\$(ls /usr/local/conda-meta/transdecoder-*.json 2>/dev/null | head -n1)" .json | sed -E 's/^transdecoder-([^-]+)-.*/\\1/')
    printf '"%s":\\n  transdecoder_predict: "%s"\\n  container: "%s"\\n' \\
      "${task.process}" "\${transdecoder_version:-unknown}" "${task.container}" > versions.yml
    """

  stub:
    """
    set -euo pipefail
    printf "mikado_stub_tx\\t0\\t42\\tmikado_stub_orf\\t0\\t+\\t0\\t42\\t0\\t1\\t42,\\t0,\\n" > mikado_prepared.fasta.transdecoder.bed
    printf ">mikado_stub_orf\\nMAIVMGR\\n" > mikado_prepared.fasta.transdecoder.pep
    cat > mikado_prepared.fasta.transdecoder.gff3 <<'EOF'
##gff-version 3
mikado_stub_tx	TransDecoder	CDS	1	21	.	+	0	ID=mikado_stub_orf;Parent=mikado_stub_tx
EOF
    printf '"%s":\n  transdecoder_predict: "stub"\n' "${task.process}" > versions.yml
    """
}
