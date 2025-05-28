process egapx {

  tag "Executing NCBI egapx gene annotation pipeline ..."
  container 'avelt/ncbi_egapx:0.3.2-alpha'
  containerOptions "--volume $egapx_paramfile_path:/egapx_param_path --volume $projectDir:$projectDir"
  publishDir "${params.output_dir}/egapx"
  cpus 5

  input:
    path(egapx_paramfile_path)
    path(egapx_paramfile_name)

  output:
    path("*")

  script:
    """
    DATE=\$(date "+%Y-%m-%d %H:%M:%S")

    CMD="python3 /opt/egapx/ui/egapx.py -o egapx_ncbi_out -e local -c /opt/egapx/egapx_config -w dfs_work -o dfs_out -lc /opt/egapx/local_cache /egapx_param_path/$egapx_paramfile_name"
    echo "[\$DATE] Executing: \$CMD"

    python3 /opt/egapx/ui/egapx.py -o egapx_ncbi_out -e local -c /opt/egapx/egapx_config -w dfs_work -o dfs_out -lc /opt/egapx/local_cache /egapx_param_path/$egapx_paramfile_name
    """
}
