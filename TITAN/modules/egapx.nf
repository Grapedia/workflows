process egapx {

  tag "Executing NCBI egapx gene annotation pipeline ..."
  container 'avelt/ncbi_egapx:0.3.2-alpha'
  containerOptions "--volume $projectDir:$projectDir"
  publishDir "${params.output_dir}/egapx"
  cpus params.egapx_cpus

  input:
    path egapx_paramfile

  output:
    path("*"), emit: egapx_results

  script:
    """
    DATE=\$(date "+%Y-%m-%d %H:%M:%S")

    CMD="python3 /opt/egapx/ui/egapx.py -e local -c /opt/egapx/work/egapx_config -w dfs_work -o \$PWD -lc /opt/egapx/local_cache ${egapx_paramfile}"
    echo "[\$DATE] Executing: \$CMD"

    python3 /opt/egapx/ui/egapx.py -e local -c /opt/egapx/work/egapx_config -w dfs_work -o \$PWD -lc /opt/egapx/local_cache ${egapx_paramfile}
    """

  stub:
    """
    mkdir -p egapx_stub
    printf "##gff-version 3\\n" > egapx_stub/egapx_annotation.gff3
    """
}
