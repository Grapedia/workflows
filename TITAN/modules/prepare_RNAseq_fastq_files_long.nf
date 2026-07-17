// We check that if the samples are of type FASTQ they exist in params.RNAseq_data_dir.
// If the samples are SRA, we'll resolve ENA FASTQ URLs and download gzipped FASTQ files.
// Once done, the RNAseq samples are ready for transcriptome assembly.

process prepare_RNAseq_fastq_files_long {
  label 'process_low'
  tag "prepare_RNAseq_fastq_files on $sample_ID"
  container params.container_python
  debug true

  input:
  tuple val(sample_ID), val(SRA_or_FASTQ), val(library_layout), path(local_reads)
  val(ena_download_timeout_seconds)
  val(ena_max_download_attempts)
  val(ena_retry_wait_seconds)
  val(ena_verify_md5)

  output:
  tuple val(sample_ID), val(SRA_or_FASTQ), val(library_layout), path("long_read_input.*"), emit: prepared_fastqs
    path "versions.yml", emit: versions



  script:
  """
  set -euo pipefail
  DATE=\$(date "+%Y-%m-%d %H:%M:%S")
  echo "[\$DATE] Preparing long RNA-seq files for $sample_ID"

  if [[ $SRA_or_FASTQ == "SRA" ]]
  then
    echo "[\$DATE] Downloading SRA long-read sample through ENA API: $sample_ID"
    python3 ${projectDir}/scripts/download_sra_fastq.py \\
      --accession $sample_ID \\
      --layout long \\
      --outdir . \\
      --timeout-seconds ${ena_download_timeout_seconds} \\
      --max-attempts ${ena_max_download_attempts} \\
      --retry-wait-seconds ${ena_retry_wait_seconds} \\
      ${ena_verify_md5 ? '--verify-md5' : '--no-verify-md5'}
  elif [[ $SRA_or_FASTQ == "FASTQ" ]]
  then
    test -s ${sample_ID}.fastq.gz
  elif [[ $SRA_or_FASTQ == "FASTA" ]]
  then
    test -s ${sample_ID}.fasta
  else
    echo "[\$DATE] ERROR: \$SRA_or_FASTQ is not equal to SRA, FASTQ or FASTA" >&2
    exit 1
  fi

  if [[ -s ${sample_ID}.fastq.gz ]]; then
    cp ${sample_ID}.fastq.gz long_read_input.fastq.gz
  elif compgen -G "${sample_ID}*.fastq.gz" > /dev/null; then
    first_fastq=\$(compgen -G "${sample_ID}*.fastq.gz" | sort | head -n 1)
    cp "\$first_fastq" long_read_input.fastq.gz
  elif [[ -s ${sample_ID}.fasta ]]; then
    cp ${sample_ID}.fasta long_read_input.fasta
  else
    echo "[\$DATE] ERROR: no prepared long-read FASTQ/FASTA found for $sample_ID" >&2
    exit 1
  fi
    printf '"%s":\n  container: "not_recorded"\n' "${task.process}" > versions.yml
  """

  stub:
  """
  rm -f ${sample_ID}*
  printf "@stub\\nACGT\\n+\\n!!!!\\n" | gzip -c > long_read_input.fastq.gz
  printf '"%s":\n  container: "not_recorded"\n' "${task.process}" > versions.yml
  """
}
