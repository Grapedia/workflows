// We check that if the samples are of type FASTQ they exist in params.RNAseq_data_dir.
// If the samples are SRA, we'll resolve ENA FASTQ URLs and download gzipped FASTQ files.
// Once done, the RNAseq samples are ready for transcriptome assembly.

process prepare_RNAseq_fastq_files_long {
  label 'process_low'
  tag "prepare_RNAseq_fastq_files on $sample_ID"
  container params.container_python
  debug true
  // ENA downloads have script-level retries and keep a process retry as an outer safety net.
  errorStrategy 'retry'

  input:
  tuple val(sample_ID), val(SRA_or_FASTQ), val(library_layout), path(local_reads)

  output:
  tuple val(sample_ID), val(SRA_or_FASTQ), val(library_layout), path("${sample_ID}*", includeInputs: true), emit : prepared_fastqs

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
      --timeout-seconds ${params.ena_download_timeout_seconds} \\
      --max-attempts ${params.ena_max_download_attempts} \\
      --retry-wait-seconds ${params.ena_retry_wait_seconds} \\
      ${params.ena_verify_md5 ? '--verify-md5' : '--no-verify-md5'}
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
  """

  stub:
  """
  rm -f ${sample_ID}*
  printf "@stub\\nACGT\\n+\\n!!!!\\n" | gzip -c > ${sample_ID}.fastq.gz
  """
}
