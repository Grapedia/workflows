process salmon_strand_inference {

  tag "Executing salmon strand inference on $cds_fasta"
  container 'quay.io/biocontainers/salmon:1.10.3--haf24da9_3'
  publishDir "$params.outdir/salmon_strand/"
  cpus 4

  input:
    tuple val(sample_ID), val(library_layout), path(reads)
    path(salmon_index)

  output:
    tuple val(sample_ID), val(library_layout), path(reads), path("${sample_ID}.strand_info"), emit : salmon_strand_info

  script:
    """
    if [[ $library_layout == "paired" ]]
    then
      salmon quant -i $salmon_index -l A -p ${task.cpus} -1 ${reads[0]} -2 ${reads[1]} -o ${sample_ID}_quant --validateMappings --skipQuant 2> ${sample_ID}.log
    elif [[ $library_layout == "single" ]]
    then
      salmon quant -i $salmon_index -l A -p ${task.cpus} -r ${reads} -o ${sample_ID}_quant --validateMappings --skipQuant 2> ${sample_ID}.log
    fi

    if grep 'Automatically detected most likely library type' ${sample_ID}.log > /dev/null
    then 
      grep 'Automatically detected most likely library type' ${sample_ID}.log | sed 's/.*type as //' > ${sample_ID}.strand_info
    else 
      echo 'U' > ${sample_ID}.strand_info
    fi
    """
}