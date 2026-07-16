nextflow.enable.dsl=2

// Include various processing modules
include { aegis_short_reads } from "../modules/aegis_short_reads"
include { aegis_long_reads } from "../modules/aegis_long_reads"
include { diamond2go } from "../modules/diamond2go"

workflow aegis {

  take:
    input_data
    run_edta
    use_long_reads

  main:

    def workflow_inputs = input_data // Create a new variable
    def edta_enabled = run_edta
    def use_long_reads_enabled = use_long_reads
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
        if (use_long_reads_enabled) {
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


    // ----------------------------------------------------------------------------------------
    //     Aegis scripts (1, 2, 3) to create the final GFF3 file from all the evidences
    // ----------------------------------------------------------------------------------------

    if (edta_enabled) {
      if (use_long_reads_enabled) {
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
    } else {
      println "Skipping the AEGIS process because EDTA was not run (run_edta = '${edta_enabled}'). To enable EDTA, and consequently AEGIS, set run_edta = true or EDTA = 'yes'."
    }

    // ----------------------------------------------------------------------------------------
    //               Diamond2GO on proteins predicted with TITAN/Aegis
    // ----------------------------------------------------------------------------------------

    if (edta_enabled) {
      if (use_long_reads_enabled) {
         diamond2go(aegis_long_reads.out.aegis_proteins_all, aegis_long_reads.out.aegis_proteins_main)
      } else {
         diamond2go(aegis_short_reads.out.aegis_proteins_all, aegis_short_reads.out.aegis_proteins_main)
      }
    } else {
      println "Skipping the Diamond2GO process because EDTA was not run (run_edta = '${edta_enabled}'). To enable EDTA, and consequently AEGIS and Diamond2GO, set run_edta = true or EDTA = 'yes'."
    }

  // Outputs
  emit:
    aegis_gff            = edta_enabled ? (use_long_reads_enabled ? aegis_long_reads.out.aegis_gff : aegis_short_reads.out.aegis_gff) : Channel.empty()
    aegis_proteins_all   = edta_enabled ? (use_long_reads_enabled ? aegis_long_reads.out.aegis_proteins_all : aegis_short_reads.out.aegis_proteins_all) : Channel.empty()
    aegis_proteins_main  = edta_enabled ? (use_long_reads_enabled ? aegis_long_reads.out.aegis_proteins_main : aegis_short_reads.out.aegis_proteins_main) : Channel.empty()
    diamond2go_output    = edta_enabled ? diamond2go.out : Channel.empty()
}
