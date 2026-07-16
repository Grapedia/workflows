nextflow.enable.dsl = 2

params.RNAseq_samplesheet = false
params.protein_samplesheet = false
params.new_assembly = false
params.previous_assembly = false
params.previous_annotations = false
params.output_dir = false
params.egapx_paramfile = false
params.egapx_version = params.egapx_version ?: "0.5.2"
params.egapx_revision = params.egapx_revision ?: "v0.5.2"
params.egapx_container = params.egapx_container ?: params.container_egapx
params.egapx_executor = params.egapx_executor ?: "docker"
params.egapx_data_version = params.egapx_data_version ?: "current_1"
params.egapx_runner_dir = params.egapx_runner_dir ?: false
params.aegis_version = params.aegis_version ?: "v0.3.25"
params.aegis_container = params.aegis_container ?: params.container_aegis

include { TITAN } from './workflows/titan'

workflow {
    TITAN()
}
