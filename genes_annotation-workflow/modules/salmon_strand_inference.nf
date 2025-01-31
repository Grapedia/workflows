process salmon_strand_inference {

  tag "Executing salmon strand inference on $sample_ID"
  container 'quay.io/biocontainers/salmon:1.10.3--haf24da9_3'
  containerOptions "--volume ${projectDir}/work:/work --volume ${projectDir}/scripts:/scripts"
  publishDir "$params.outdir/salmon_strand/"
  cpus 4

  input:
    tuple val(sample_ID), val(library_layout), path(reads)
    path(salmon_index)

  output:
    tuple val(sample_ID), val(library_layout), path(reads), path("${sample_ID}.strand_info.classified")

  script:
    """
    DATE=\$(date "+%Y-%m-%d %H:%M:%S")
    echo "[\$DATE] Running salmon strand inference on $sample_ID" >> ${params.logfile} 2>&1

    index_path=\$(/scripts/retrieve_path.sh $salmon_index)

    if [[ $library_layout == "paired" ]]
    then
      CMD="salmon quant -i \${index_path} -l A -p ${task.cpus} -1 ${reads[0]} -2 ${reads[1]} -o ${sample_ID}_quant --validateMappings --skipQuant 2> ${sample_ID}.log"
      echo "[\$DATE] Executing: \$CMD" >> ${params.logfile} 2>&1
      salmon quant -i \${index_path} -l A -p ${task.cpus} -1 ${reads[0]} -2 ${reads[1]} -o ${sample_ID}_quant --validateMappings --skipQuant 2> ${sample_ID}.log
    elif [[ $library_layout == "single" ]]
    then
      CMD="salmon quant -i \${index_path} -l A -p ${task.cpus} -r ${reads} -o ${sample_ID}_quant --validateMappings --skipQuant 2> ${sample_ID}.log"
      echo "[\$DATE] Executing: \$CMD" >> ${params.logfile} 2>&1
      salmon quant -i \${index_path} -l A -p ${task.cpus} -r ${reads} -o ${sample_ID}_quant --validateMappings --skipQuant 2> ${sample_ID}.log
    fi

    if grep 'Automatically detected most likely library type' ${sample_ID}.log > /dev/null
    then 
      grep 'Automatically detected most likely library type' ${sample_ID}.log | sed 's/.*type as //' > ${sample_ID}.strand_info
    else 
      echo 'U' > ${sample_ID}.strand_info
    fi

    strand_type=\$(cat "${sample_ID}.strand_info")

    if [[ "\$strand_type" == "IU" || "\$strand_type" == "U" ]]; then
        strand_info="unstranded"
        echo "[\$DATE] $sample_ID is classified as unstranded" >> ${params.logfile} 2>&1
    elif [[ "\$strand_type" == "ISR" || "\$strand_type" == "FR" ]]; then
        strand_info="stranded_forward"
        echo "[\$DATE] $sample_ID is classified as stranded_forward" >> ${params.logfile} 2>&1
    elif [[ "\$strand_type" == "ISF" || "\$strand_type" == "RF" ]]; then
        strand_info="stranded_reverse"
        echo "[\$DATE] $sample_ID is classified as stranded_reverse" >> ${params.logfile} 2>&1
    else
        strand_info="unstranded"
        echo "[\$DATE] $sample_ID is classified as unstranded" >> ${params.logfile} 2>&1
    fi

    echo "\$strand_info" > "${sample_ID}.strand_info.classified"
    """
}