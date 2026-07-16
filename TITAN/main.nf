nextflow.enable.dsl = 2

include { TITAN } from './workflows/titan'

workflow {
    TITAN()
}
