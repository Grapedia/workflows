// We check that if the samples are of type FASTQ they exist in params.RNAseq_data_dir.
// If the samples are SRA, we'll resolve ENA FASTQ URLs and download gzipped FASTQ files.
// Once done, the RNAseq samples are ready for transcriptome assembly.

process prepare_RNAseq_fastq_files_short {
  label 'process_low'
  tag "prepare_RNAseq_fastq_files on $sample_ID"
  container params.container_python
  debug true
  // ENA downloads have script-level retries and keep a process retry as an outer safety net.
  errorStrategy 'retry'

  input:
  tuple val(sample_ID), val(SRA_or_FASTQ), val(library_layout), path(local_reads)
  val(ena_download_timeout_seconds)
  val(ena_max_download_attempts)
  val(ena_retry_wait_seconds)
  val(ena_verify_md5)

  output:
  tuple val(sample_ID), val(SRA_or_FASTQ), val(library_layout), path("${sample_ID}*.fastq.gz", includeInputs: true), emit: prepared_fastqs
    path "versions.yml", emit: versions



  script:
  """
  set -euo pipefail
  DATE=\$(date "+%Y-%m-%d %H:%M:%S")
  echo "[\$DATE] Preparing short RNA-seq files for $sample_ID"

  if [[ $SRA_or_FASTQ == "SRA" ]]
  then
    if [[ $library_layout == "paired" ]]
    then
      echo "[\$DATE] Downloading SRA paired-end sample through ENA API: $sample_ID"
      python3 ${projectDir}/scripts/download_sra_fastq.py \\
        --accession $sample_ID \\
        --layout paired \\
        --outdir . \\
        --timeout-seconds ${ena_download_timeout_seconds} \\
        --max-attempts ${ena_max_download_attempts} \\
        --retry-wait-seconds ${ena_retry_wait_seconds} \\
        ${ena_verify_md5 ? '--verify-md5' : '--no-verify-md5'}
    elif [[ $library_layout == "single" ]]
    then
      echo "[\$DATE] Downloading SRA single-end sample through ENA API: $sample_ID"
      python3 ${projectDir}/scripts/download_sra_fastq.py \\
        --accession $sample_ID \\
        --layout single \\
        --outdir . \\
        --timeout-seconds ${ena_download_timeout_seconds} \\
        --max-attempts ${ena_max_download_attempts} \\
        --retry-wait-seconds ${ena_retry_wait_seconds} \\
        ${ena_verify_md5 ? '--verify-md5' : '--no-verify-md5'}
    else
      echo "[\$DATE] ERROR: $library_layout is not equal to paired or single" >&2
      exit 1
    fi
  elif [[ $SRA_or_FASTQ == "FASTQ" ]]
  then
    if [[ $library_layout == "paired" ]]; then
      test -s ${sample_ID}_1.fastq.gz
      test -s ${sample_ID}_2.fastq.gz
    elif [[ $library_layout == "single" ]]; then
      test -s ${sample_ID}.fastq.gz
    else
      echo "[\$DATE] ERROR: $library_layout is not equal to paired or single" >&2
      exit 1
    fi
  else
    echo "[\$DATE] ERROR: \$SRA_or_FASTQ is not equal to SRA or FASTQ" >&2
    exit 1
  fi
    printf '"%s":\n  container: "not_recorded"\n' "${task.process}" > versions.yml
  """

  stub:
  """
  rm -f ${sample_ID}*.fastq.gz
  if [[ $library_layout == "paired" ]]
  then
    printf "@stub/1\\nACGT\\n+\\n!!!!\\n" | gzip -c > ${sample_ID}_1.fastq.gz
    printf "@stub/2\\nTGCA\\n+\\n!!!!\\n" | gzip -c > ${sample_ID}_2.fastq.gz
  else
    printf "@stub\\nACGT\\n+\\n!!!!\\n" | gzip -c > ${sample_ID}.fastq.gz
  fi
    printf '"%s":\n  container: "not_recorded"\n' "${task.process}" > versions.yml
  """
}
