// 1. Generate reference genome indices using STAR
process star_genome_indices {
  label 'process_index'

  tag "STAR genomeGenerate on $genome"
  container params.container_star
  publishDir "${params.output_dir}/intermediate_files/evidence_data/star_databases/", mode: "copy", enabled: params.publish_intermediates
  input:
    path(genome_fasta)
    val(genome)
    val(star_genome_sa_index_nbases)
    val(star_sjdb_gtf_file)

  output:
    path("${genome}_index"), emit: index
    path "versions.yml", emit: versions



  script:
    """
    set -euo pipefail
    DATE=\$(date "+%Y-%m-%d %H:%M:%S")
    echo "[\$DATE] Running STAR index creation on $genome"

    GENOME_SA_INDEX_NBASES="${star_genome_sa_index_nbases}"
    if [[ -z "\$GENOME_SA_INDEX_NBASES" || "\$GENOME_SA_INDEX_NBASES" == "false" ]]; then
      genome_bases=\$(awk '/^>/ { next } { gsub(/[ \\t\\r]/, ""); n += length(\$0) } END { print n + 0 }' ${genome_fasta})
      if [[ "\$genome_bases" -le 0 ]]; then
        echo "[\$DATE] ERROR: cannot derive genomeSAindexNbases from empty FASTA ${genome_fasta}" >&2
        exit 1
      fi
      GENOME_SA_INDEX_NBASES=\$(awk -v n="\$genome_bases" 'BEGIN { v = int(log(n) / log(2) / 2 - 1); if (v < 1) v = 1; if (v > 14) v = 14; print v }')
    elif ! [[ "\$GENOME_SA_INDEX_NBASES" =~ ^[0-9]+$ && "\$GENOME_SA_INDEX_NBASES" -ge 1 && "\$GENOME_SA_INDEX_NBASES" -le 14 ]]; then
      echo "[\$DATE] ERROR: STAR_genomeSAindexNbases must be false or an integer between 1 and 14, got '\$GENOME_SA_INDEX_NBASES'" >&2
      exit 1
    fi

    SJDB_GTF="${star_sjdb_gtf_file}"
    sjdb_args=()
    if [[ -n "\$SJDB_GTF" && "\$SJDB_GTF" != "false" ]]; then
      if [[ ! -s "\$SJDB_GTF" ]]; then
        echo "[\$DATE] ERROR: STAR_sjdbGTFfile does not exist or is empty: \$SJDB_GTF" >&2
        exit 1
      fi
      sjdb_args=(--sjdbGTFfile "\$SJDB_GTF")
    fi

    mkdir -p ${genome}_index
    CMD=(STAR --runThreadN ${task.cpus} --runMode genomeGenerate --genomeDir ${genome}_index --genomeFastaFiles ${genome_fasta} --genomeSAindexNbases "\$GENOME_SA_INDEX_NBASES" "\${sjdb_args[@]}")
    echo "[\$DATE] Executing: \${CMD[*]}"
    "\${CMD[@]}"

    STAR --version 2>&1 | sed 's/^/  star: "/; s/\$/"/' | {
      printf '"%s":\\n  genomeSAindexNbases: "%s"\\n' "${task.process}" "\$GENOME_SA_INDEX_NBASES"
      cat
    } > versions.yml
    """

  stub:
    """
    set -euo pipefail

    mkdir -p ${genome}_index
    printf "stub STAR index\\n" > ${genome}_index/Genome
    printf '"%s":\n  star: "stub"\n  genomeSAindexNbases: "stub"\n' "${task.process}" > versions.yml
    """
}
