process salmon_index {

  tag "Executing salmon indexing on $cds_fasta"
  container 'quay.io/biocontainers/salmon:1.10.3--haf24da9_3'
  containerOptions "--volume ${projectDir}/work:/work --volume ${projectDir}/scripts:/scripts"
  publishDir "$params.outdir/salmon_index/"
  cpus 4

  input:
    val(cds_fasta)

  output:
    path("salmon_index")

  script:
    """
    DATE=\$(date "+%Y-%m-%d %H:%M:%S")
    echo "[\$DATE] Running salmon indexing on ${cds_fasta}" >> ${params.logfile} 2>&1
    index_path=\$(/scripts/retrieve_path.sh ${cds_fasta})
    CMD="salmon index -t \${index_path} -i salmon_index"
    echo "[\$DATE] Executing: \$CMD" >> ${params.logfile} 2>&1
    salmon index -t \${index_path} -i salmon_index
    """
}
