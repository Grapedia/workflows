process salmon_strand_inference {

  label 'process_low'
  tag "Executing salmon strand inference on $sample_ID"
  container params.container_salmon
  // Avoid republishing the trimmed FASTQs re-emitted for downstream alignment.
  publishDir "${params.output_dir}/intermediate_files/salmon_strand/", mode: "copy", enabled: params.publish_intermediates, saveAs: { filename ->
    filename.endsWith('.trimmed.fastq.gz') ? null : filename
  }
  input:
    tuple val(sample_ID), val(library_layout), path(read_1), path(read_2)
    path(salmon_index)

  output:
    tuple val(sample_ID), val(library_layout), path(read_1), path(read_2), env('STRAND_INFO'), emit: strand_inference_result
    path("${sample_ID}.strand_info.classified"), emit: strand_classification
    path("${sample_ID}.log"), emit: salmon_log
    path "versions.yml", emit: versions



  script:
    """
    set -euo pipefail
    DATE=\$(date "+%Y-%m-%d %H:%M:%S")
    echo "[\$DATE] Running salmon strand inference on $sample_ID"

    if [[ "${library_layout}" == "paired" ]]; then
      CMD="salmon quant -i ${salmon_index} -l A -p ${task.cpus} -1 ${read_1} -2 ${read_2} -o ${sample_ID}_quant --validateMappings --skipQuant 2> ${sample_ID}.log"
      echo "[\$DATE] Executing: \$CMD"
      salmon quant -i ${salmon_index} -l A -p ${task.cpus} -1 ${read_1} -2 ${read_2} -o ${sample_ID}_quant --validateMappings --skipQuant 2> ${sample_ID}.log
    elif [[ "${library_layout}" == "single" ]]; then
      CMD="salmon quant -i ${salmon_index} -l A -p ${task.cpus} -r ${read_1} -o ${sample_ID}_quant --validateMappings --skipQuant 2> ${sample_ID}.log"
      echo "[\$DATE] Executing: \$CMD"
      salmon quant -i ${salmon_index} -l A -p ${task.cpus} -r ${read_1} -o ${sample_ID}_quant --validateMappings --skipQuant 2> ${sample_ID}.log
    else
      echo "ERROR: Unsupported library_layout '${library_layout}' for sample '${sample_ID}'. Expected 'single' or 'paired'." >&2
      exit 1
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
        echo "[\$DATE] $sample_ID is classified as unstranded"
    elif [[ "\$strand_type" == "ISF" || "\$strand_type" == "SF" || "\$strand_type" == "FR" ]]; then
        strand_info="stranded_forward"
        echo "[\$DATE] $sample_ID is classified as stranded_forward"
    elif [[ "\$strand_type" == "ISR" || "\$strand_type" == "SR" || "\$strand_type" == "RF" ]]; then
        strand_info="stranded_reverse"
        echo "[\$DATE] $sample_ID is classified as stranded_reverse"
    else
        strand_info="unstranded"
        echo "[\$DATE] $sample_ID is classified as unstranded"
    fi

    echo "\$strand_info" > "${sample_ID}.strand_info.classified"
    export STRAND_INFO="\$strand_info"
    salmon --version | sed 's/^/  salmon: "/; s/\$/"/' | {
      printf '"%s":\\n' "${task.process}"
      cat
    } > versions.yml
    """

  stub:
    """
    set -euo pipefail

    export STRAND_INFO="stranded_forward"
    echo "\$STRAND_INFO" > ${sample_ID}.strand_info.classified
    printf "stub salmon log\\n" > ${sample_ID}.log
    printf '"%s":\n  salmon: "stub"\n' "${task.process}" > versions.yml
    """
}
