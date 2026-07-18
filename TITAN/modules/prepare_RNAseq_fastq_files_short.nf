// Prepare short-read RNA-seq inputs from local FASTQ files or ENA-resolved SRA accessions.

process prepare_RNAseq_fastq_files_short {
  label 'process_low'
  tag "prepare_RNAseq_fastq_files on $sample_ID"
  container params.container_python

  input:
  tuple val(sample_ID), val(SRA_or_FASTQ), val(library_layout), path(local_reads)
  path(download_sra_fastq_script)
  val(ena_download_timeout_seconds)
  val(ena_max_download_attempts)
  val(ena_retry_wait_seconds)
  val(ena_verify_md5)

  output:
  tuple val(sample_ID), val(SRA_or_FASTQ), val(library_layout), path("prepared_1.fastq.gz"), path("prepared_2.fastq.gz"), emit: prepared_fastqs
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
      python3 ${download_sra_fastq_script} \\
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
      python3 ${download_sra_fastq_script} \\
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

  if [[ $library_layout == "paired" ]]; then
    cp ${sample_ID}_1.fastq.gz prepared_1.fastq.gz
    cp ${sample_ID}_2.fastq.gz prepared_2.fastq.gz
  else
    cp ${sample_ID}.fastq.gz prepared_1.fastq.gz
    printf "" | gzip -c > prepared_2.fastq.gz
  fi

  python3 --version 2>&1 | sed 's/^/  python: "/; s/\$/"/' | {
    printf '"%s":\\n' "${task.process}"
    cat
  } > versions.yml
  """

  stub:
  """
  set -euo pipefail

  if [[ $library_layout == "paired" ]]
  then
    printf "@stub/1\\nACGT\\n+\\n!!!!\\n" | gzip -c > prepared_1.fastq.gz
    printf "@stub/2\\nTGCA\\n+\\n!!!!\\n" | gzip -c > prepared_2.fastq.gz
  else
    printf "@stub\\nACGT\\n+\\n!!!!\\n" | gzip -c > prepared_1.fastq.gz
    printf "" | gzip -c > prepared_2.fastq.gz
  fi
  printf '"%s":\n  python: "stub"\n' "${task.process}" > versions.yml
  """
}
