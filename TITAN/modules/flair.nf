process flair_isoforms {
  label 'process_transcriptome'
  tag "FLAIR isoforms on ${sample_ID}"
  container params.container_flair
  publishDir "${params.output_dir}/additional_annotations/flair/samples", mode: 'copy', saveAs: { filename ->
    if (filename.endsWith('.flair.isoforms.gtf') || filename.endsWith('.flair.isoforms.fa') || filename.endsWith('.flair.log') || filename == 'versions.yml') {
      return filename
    }
    return null
  }

  input:
    tuple val(sample_ID), val(SRA_or_FASTQ), val(library_layout), val(read_format), path(reads_fastq), path(reads_fasta)
    path(genome)
    path(liftoff_gtf)

  output:
    tuple val(sample_ID), path("${sample_ID}.flair.isoforms.gtf"), path("${sample_ID}.flair.isoforms.fa"), emit: isoforms
    path "${sample_ID}.flair.log", emit: log
    path "versions.yml", emit: versions

  script:
    """
    set -euo pipefail
    if [[ "${params.run_flair}" != "true" ]]; then
      : > "${sample_ID}.flair.isoforms.gtf"
      : > "${sample_ID}.flair.isoforms.fa"
      printf "FLAIR skipped; set --run_flair true to enable.\\n" > "${sample_ID}.flair.log"
      printf '"%s":\n  flair: "skipped"\n  container: "%s"\n' "${task.process}" "${task.container}" > versions.yml
      exit 0
    fi

    if [[ "${read_format}" == "fastq" ]]; then
      reads="${reads_fastq}"
    elif [[ "${read_format}" == "fasta" ]]; then
      reads="${reads_fasta}"
    else
      echo "ERROR: read_format must be fastq or fasta, got ${read_format}" >&2
      exit 1
    fi

    # flair correct has no genome (-g) argument in this pinned version, and -f/--gtf
    # strictly requires GTF (liftoff_gtf is pre-converted from GFF3 upstream; see
    # modules/agat_convert_gff3_to_gtf.nf).
    {
      flair align -g "${genome}" -r "\${reads}" -o "${sample_ID}.flair" --threads ${task.cpus}
      flair correct -q "${sample_ID}.flair.bed" -f "${liftoff_gtf}" -o "${sample_ID}.flair.corrected" --threads ${task.cpus}
      flair collapse -g "${genome}" -r "\${reads}" -q "${sample_ID}.flair.corrected_all_corrected.bed" -o "${sample_ID}.flair" --threads ${task.cpus}
    } 2>&1 | tee "${sample_ID}.flair.log"

    test -s "${sample_ID}.flair.isoforms.gtf"
    test -s "${sample_ID}.flair.isoforms.fa"
    flair --version 2>&1 | sed 's/^/  flair: "/; s/\$/"/' | {
      printf '"%s":\\n' "${task.process}"
      cat
      printf '  container: "%s"\\n' "${task.container}"
    } > versions.yml
    """

  stub:
    """
    set -euo pipefail
    seqid=\$(awk '/^>/ {print substr(\$1, 2); exit}' ${genome})
    end=\$(awk 'BEGIN{n=0} !/^>/ {n += length(\$0)} END{print n < 30 ? n : 30}' ${genome})
    printf "%s\\tFLAIR\\ttranscript\\t1\\t%s\\t.\\t+\\t.\\tgene_id \\"%s_flair_gene\\"; transcript_id \\"%s_flair_tx\\";\\n" "\$seqid" "\$end" "${sample_ID}" "${sample_ID}" > "${sample_ID}.flair.isoforms.gtf"
    printf ">%s_flair_tx\\nATGGCCATTGTAATGGGCCGCTGAAAG\\n" "${sample_ID}" > "${sample_ID}.flair.isoforms.fa"
    printf "FLAIR stub\\n" > "${sample_ID}.flair.log"
    printf '"%s":\n  flair: "stub"\n' "${task.process}" > versions.yml
    """
}

process flair_merge_isoforms {
  label 'process_merge'
  tag "Merge FLAIR isoform GTFs"
  container params.container_python
  publishDir "${params.output_dir}/additional_annotations/flair", mode: 'copy', saveAs: { filename ->
    if (filename in ['flair_isoforms.gtf', 'flair_isoforms.fa', 'versions.yml']) {
      return filename
    }
    return null
  }

  input:
    path(flair_gtfs, stageAs: "flair_gtfs/*")
    path(flair_fastas, stageAs: "flair_fastas/*")

  output:
    path "flair_isoforms.gtf", emit: gtf
    path "flair_isoforms.fa", emit: fasta
    path "versions.yml", emit: versions

  script:
    """
    set -euo pipefail
    : > flair_isoforms.gtf
    : > flair_isoforms.fa
    find flair_gtfs -type f -name '*.gtf' -size +0c -print0 | sort -z | xargs -0 -r cat >> flair_isoforms.gtf
    find flair_fastas -type f -name '*.fa' -size +0c -print0 | sort -z | xargs -0 -r cat >> flair_isoforms.fa
    printf '"%s":\n  flair_merge_isoforms: "python-stdlib"\n  container: "%s"\n' "${task.process}" "${task.container}" > versions.yml
    """

  stub:
    """
    set -euo pipefail
    cat flair_gtfs/*.gtf > flair_isoforms.gtf 2>/dev/null || : > flair_isoforms.gtf
    cat flair_fastas/*.fa > flair_isoforms.fa 2>/dev/null || : > flair_isoforms.fa
    printf '"%s":\n  flair_merge_isoforms: "stub"\n' "${task.process}" > versions.yml
    """
}

process flair_empty_evidence {
  label 'process_low'
  tag "Create empty FLAIR evidence sentinels"
  container params.container_python

  input:
    path(empty_gtf)

  output:
    path "flair_isoforms.gtf", emit: gtf
    path "flair_isoforms.fa", emit: fasta
    path "versions.yml", emit: versions

  script:
    """
    set -euo pipefail
    cp "${empty_gtf}" flair_isoforms.gtf
    : > flair_isoforms.fa
    printf '"%s":\n  flair: "empty"\n  container: "%s"\n' "${task.process}" "${task.container}" > versions.yml
    """

  stub:
    """
    set -euo pipefail
    : > flair_isoforms.gtf
    : > flair_isoforms.fa
    printf '"%s":\n  flair: "empty_stub"\n' "${task.process}" > versions.yml
    """
}
