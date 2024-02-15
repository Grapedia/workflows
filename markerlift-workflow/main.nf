nextflow.enable.dsl=2

params.genome_old = "$projectDir/genome_old.fasta"
params.genome_new = "$projectDir/genome_new.fasta"
params.coords_old = "$projectDir/coords_old.vcf"
params.coords_new = "$projectDir/coords_new.vcf"
params.window_lenght = 300

include { SNPLIFT } from "./modules/snplift"
include { SNPLIFT as SNPLIFT_START } from "./modules/snplift"
include { SNPLIFT as SNPLIFT_END } from "./modules/snplift"
include { SPLIT_COORDS } from "./modules/split_coords" 
include { CHECK_COORDS } from "./modules/check_coords" 
include { MERGE_RESULT } from "./modules/merge_result"
include { CREATE_CONFIG_ONLY_START } from "./modules/only_start"
include { CREATE_CONFIG_START; CREATE_CONFIG_END } from "./modules/start_end"

log.info """
     Markerlift is heavily based on
     SNPlift: Lift over SNP positions to match a new reference genome
     Source: https://github.com/enormandeau/snplift

     Coords input file is a tab delimited file with 4 columns
     <chr> <start> <end> <name>

     if <start> is equal to <end> it just run as SNPlift
     Otherwise it runs SNPlift on <start> then it runs 
     SNPlift on <end> and then merges the results.
     =====================================================
     ----- INPUT
     old genome  : ${params.genome_old}
     new genome  : ${params.genome_new}
     coords file : ${params.coords_old}
     window length: ${params.window_lenght}
     ---- OUTPUT
     output file : ${params.coords_new}
     """
     .stripIndent()

workflow ONLY_START {
    take:
      first_ch
    main:
      config_file = CREATE_CONFIG_ONLY_START(first_ch,
        params.genome_old,
        params.genome_new,
        params.coords_old,
        params.coords_new,
        params.window_lenght
      )
      SNPLIFT(config_file, "$projectDir/config_only_start.sh", params.coords_new)
}

workflow START_END {
    take:
      second_ch
    main:
      split = SPLIT_COORDS("$projectDir/$params.coords_old")
      config_start = CREATE_CONFIG_START(second_ch,
        params.genome_old,
        params.genome_new,
        split.start_coords,
        "START_$params.coords_new",
        params.window_lenght
      )
      config_end = CREATE_CONFIG_END(second_ch,
        params.genome_old,
        params.genome_new,
        split.end_coords,
        "END_$params.coords_new",
        params.window_lenght
      )
      config_1 = SNPLIFT_START(config_start, "$projectDir/config_start.sh", split.start_coords)
      config_2 = SNPLIFT_END(config_1, "$projectDir/config_end.sh", split.end_coords)
      MERGE_RESULT(config_2, "$projectDir/START_$params.coords_new", "$projectDir/END_$params.coords_new", "$projectDir/$params.coords_new")
}

workflow{
    coords_type = CHECK_COORDS("$projectDir/$params.coords_old").toInteger()
    
    coords_type.branch {
        only_start: it == 1
        start_end: it == 0
    }.set { result }

    ONLY_START(result.only_start)
    START_END(result.start_end)
}