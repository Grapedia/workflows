nextflow.enable.dsl=2

// Include various processing modules
include { aegis_short_reads } from "../modules/aegis_short_reads"
include { aegis_long_reads } from "../modules/aegis_long_reads"
include { diamond2go } from "../modules/diamond2go"

workflow aegis {

  take:
    input_data
    has_long_reads

  main:

    def workflow_inputs = input_data // Create a new variable
    def has_long_reads_enabled = has_long_reads
    // println "Files sent to Aegis process : ${workflow_inputs}"

    def masked_genome = []
    def braker3_prediction_augustus = []
    def Stringtie_merging_long_reads_default_args = []
    def Stringtie_merging_long_reads_alt_args = []
    def braker3_prediction_genemark = []
    def liftoff_annotation = []
    def Stringtie_merging_short_reads_STAR_default_args_stranded = []
    def Stringtie_merging_short_reads_STAR_alt_args_stranded = []
    def Stringtie_merging_short_reads_STAR_default_args_unstranded = []
    def Stringtie_merging_short_reads_STAR_alt_args_unstranded = []
    def gffcompare_star_psiclass_stranded = []
    def gffcompare_star_psiclass_unstranded = []

    workflow_inputs.each { entry ->
        def key = entry[0]
        def filepath = entry[1]
        if (key == "masked_genome.masked_genome") {
            masked_genome << filepath
            println "Masked genome file: ${masked_genome}"
        }
        if (key == "braker3_results.augustus_gff") {
            braker3_prediction_augustus << filepath
            println "braker3_prediction_augustus file: ${braker3_prediction_augustus}"
        }
        if (key == "braker3_results.genemark_gtf") {
            braker3_prediction_genemark << filepath
            println "braker3_prediction_genemark file: ${braker3_prediction_genemark}"
        }
        if (key == "previous_annotations.liftoff_previous_annotations") {
            liftoff_annotation << filepath
            println "liftoff_annotation file: ${liftoff_annotation}"
        }
        if (key == "merged_star_stringtie.default_args_stranded") {
            Stringtie_merging_short_reads_STAR_default_args_stranded << filepath
            println "Stringtie_merging_short_reads_STAR_default_args_stranded file: ${Stringtie_merging_short_reads_STAR_default_args_stranded}"
        }
        if (key == "merged_star_stringtie.alt_args_stranded") {
            Stringtie_merging_short_reads_STAR_alt_args_stranded << filepath
            println "Stringtie_merging_short_reads_STAR_alt_args_stranded file: ${Stringtie_merging_short_reads_STAR_alt_args_stranded}"
        }
        if (key == "merged_star_stringtie.default_args_unstranded") {
            if (!filepath.name.equals("dev_null1")) {
                Stringtie_merging_short_reads_STAR_default_args_unstranded << filepath
                println "Stringtie_merging_short_reads_STAR_default_args_unstranded file: ${Stringtie_merging_short_reads_STAR_default_args_unstranded}"
            } else {
                Stringtie_merging_short_reads_STAR_default_args_unstranded << filepath
                println "Stringtie_merging_short_reads_STAR_default_args_unstranded : no unstranded data available, pass ..."
            }
        }
        if (key == "merged_star_stringtie.alt_args_unstranded") {
            if (!filepath.name.equals("dev_null2")) {
                Stringtie_merging_short_reads_STAR_alt_args_unstranded << filepath
                println "Stringtie_merging_short_reads_STAR_alt_args_unstranded file: ${Stringtie_merging_short_reads_STAR_alt_args_unstranded}"
            } else {
                Stringtie_merging_short_reads_STAR_alt_args_unstranded << filepath
                println "Stringtie_merging_short_reads_STAR_alt_args_unstranded : no unstranded data available, pass ..."
            }
        }
        if (key == "gffcompare_out.star_psiclass_stranded") {
            gffcompare_star_psiclass_stranded << filepath
            println "gffcompare_star_psiclass_stranded file: ${gffcompare_star_psiclass_stranded}"
        }
        if (key == "gffcompare_out.star_psiclass_unstranded") {
            if (!filepath.name.equals("dev_null3")) {
                gffcompare_star_psiclass_unstranded << filepath
                println "gffcompare_star_psiclass_unstranded file: ${gffcompare_star_psiclass_unstranded}"
            } else {
                gffcompare_star_psiclass_unstranded << filepath
                println "gffcompare_star_psiclass_unstranded : no unstranded data available, pass ..."
            }
        }
        if (has_long_reads_enabled) {
          if (key == "merged_long_reads.default_args_gff") {
            Stringtie_merging_long_reads_default_args << filepath
            println "Stringtie_merging_long_reads_default_args file: ${Stringtie_merging_long_reads_default_args}"
          }
          if (key == "merged_long_reads.alt_args_gff") {
            Stringtie_merging_long_reads_alt_args << filepath
            println "Stringtie_merging_long_reads_alt_args file: ${Stringtie_merging_long_reads_alt_args}"
          }
        }
    }


    def missing_required_inputs = []
    [
        "masked_genome.masked_genome": masked_genome,
        "braker3_results.augustus_gff": braker3_prediction_augustus,
        "braker3_results.genemark_gtf": braker3_prediction_genemark,
        "previous_annotations.liftoff_previous_annotations": liftoff_annotation,
        "merged_star_stringtie.default_args_stranded": Stringtie_merging_short_reads_STAR_default_args_stranded,
        "merged_star_stringtie.alt_args_stranded": Stringtie_merging_short_reads_STAR_alt_args_stranded,
        "gffcompare_out.star_psiclass_stranded": gffcompare_star_psiclass_stranded,
        "merged_star_stringtie.default_args_unstranded": Stringtie_merging_short_reads_STAR_default_args_unstranded,
        "merged_star_stringtie.alt_args_unstranded": Stringtie_merging_short_reads_STAR_alt_args_unstranded,
        "gffcompare_out.star_psiclass_unstranded": gffcompare_star_psiclass_unstranded
    ].each { key, values ->
        if (!values) {
            missing_required_inputs << key
        }
    }
    if (has_long_reads_enabled) {
        [
            "merged_long_reads.default_args_gff": Stringtie_merging_long_reads_default_args,
            "merged_long_reads.alt_args_gff": Stringtie_merging_long_reads_alt_args
        ].each { key, values ->
            if (!values) {
                missing_required_inputs << key
            }
        }
    }
    if (missing_required_inputs) {
        error "Missing required Aegis input(s): ${missing_required_inputs.join(', ')}"
    }

    // ----------------------------------------------------------------------------------------
    //     Aegis scripts (1, 2, 3) to create the final GFF3 file from all the evidences
    // ----------------------------------------------------------------------------------------

    if (has_long_reads_enabled) {
      aegis_long_reads(
        file(params.new_assembly).getParent(),
        file(params.new_assembly).getName(),
        file(params.protein_samplesheet).getParent(),
        file(params.protein_samplesheet).getName(),
        masked_genome[0],
        braker3_prediction_augustus[0],
        braker3_prediction_genemark[0],
        liftoff_annotation[0],
        Stringtie_merging_long_reads_default_args[0],
        Stringtie_merging_long_reads_alt_args[0],
        Stringtie_merging_short_reads_STAR_default_args_stranded[0],
        Stringtie_merging_short_reads_STAR_alt_args_stranded[0],
        gffcompare_star_psiclass_stranded[0],
        gffcompare_star_psiclass_unstranded[0],
        Stringtie_merging_short_reads_STAR_default_args_unstranded[0],
        Stringtie_merging_short_reads_STAR_alt_args_unstranded[0]
        )
    } else {
      aegis_short_reads(
        file(params.new_assembly).getParent(),
        file(params.new_assembly).getName(),
        file(params.protein_samplesheet).getParent(),
        file(params.protein_samplesheet).getName(),
        masked_genome[0],
        braker3_prediction_augustus[0],
        braker3_prediction_genemark[0],
        liftoff_annotation[0],
        Stringtie_merging_short_reads_STAR_default_args_stranded[0],
        Stringtie_merging_short_reads_STAR_alt_args_stranded[0],
        gffcompare_star_psiclass_stranded[0],
        gffcompare_star_psiclass_unstranded[0],
        Stringtie_merging_short_reads_STAR_default_args_unstranded[0],
        Stringtie_merging_short_reads_STAR_alt_args_unstranded[0]
        )
    }

    // ----------------------------------------------------------------------------------------
    //               Diamond2GO on proteins predicted with TITAN/Aegis
    // ----------------------------------------------------------------------------------------

    if (has_long_reads_enabled) {
       diamond2go(aegis_long_reads.out.aegis_proteins_all, aegis_long_reads.out.aegis_proteins_main)
    } else {
       diamond2go(aegis_short_reads.out.aegis_proteins_all, aegis_short_reads.out.aegis_proteins_main)
    }

  // Outputs
  emit:
    aegis_gff            = has_long_reads_enabled ? aegis_long_reads.out.aegis_gff : aegis_short_reads.out.aegis_gff
    aegis_proteins_all   = has_long_reads_enabled ? aegis_long_reads.out.aegis_proteins_all : aegis_short_reads.out.aegis_proteins_all
    aegis_proteins_main  = has_long_reads_enabled ? aegis_long_reads.out.aegis_proteins_main : aegis_short_reads.out.aegis_proteins_main
    diamond2go_output    = diamond2go.out
}
