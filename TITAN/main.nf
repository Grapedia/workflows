nextflow.enable.dsl = 2

params.RNAseq_samplesheet = false
params.protein_samplesheet = false
params.new_assembly = false
params.previous_assembly = false
params.previous_annotations = false
params.output_dir = false
params.egapx_paramfile = false

params.workflow = params.workflow ?: "generate_evidence_data"

include { TITAN } from './workflows/titan'

workflow {
    TITAN()
}
