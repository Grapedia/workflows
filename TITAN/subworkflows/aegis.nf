nextflow.enable.dsl=2

// Include various processing modules
include { aegis_short_reads } from "../modules/aegis_short_reads"
include { aegis_long_reads } from "../modules/aegis_long_reads"
include { diamond2go } from "../modules/diamond2go"

workflow aegis {

  take:
    input_data

  main:

    def workflow_inputs = input_data // Create a new variable
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

    for (entry in workflow_inputs) {
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
            Stringtie_merging_short_reads_STAR_default_args_unstranded << filepath
            println "Stringtie_merging_short_reads_STAR_default_args_unstranded file: ${Stringtie_merging_short_reads_STAR_default_args_unstranded}"
        }
        if (key == "merged_star_stringtie.alt_args_unstranded") {
            Stringtie_merging_short_reads_STAR_alt_args_unstranded << filepath
            println "Stringtie_merging_short_reads_STAR_alt_args_unstranded file: ${Stringtie_merging_short_reads_STAR_alt_args_unstranded}"
        }
        if (key == "gffcompare_out.star_psiclass_stranded") {
            gffcompare_star_psiclass_stranded << filepath
            println "gffcompare_star_psiclass_stranded file: ${gffcompare_star_psiclass_stranded}"
        }
        if (key == "gffcompare_out.star_psiclass_unstranded") {
            gffcompare_star_psiclass_unstranded << filepath
            println "gffcompare_star_psiclass_unstranded file: ${gffcompare_star_psiclass_unstranded}"
        }
        if (params.use_long_reads) {
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

    if (params.EDTA == 'yes') {
      if (params.use_long_reads) {
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
      println "Skipping the AEGIS process because EDTA was not run (params.EDTA = '${params.EDTA}'). To enable EDTA, and consequently AEGIS, set EDTA = 'yes' in the nextflow.config file."
    }

    // ----------------------------------------------------------------------------------------
    //               Diamond2GO on proteins predicted with TITAN/Aegis
    // ----------------------------------------------------------------------------------------
    
    if (params.EDTA == 'yes') {
      if (params.use_long_reads) {
         diamond2go(aegis_long_reads.out.aegis_proteins_main, aegis_long_reads.out.aegis_proteins_all)
      } else {
         diamond2go(aegis_short_reads.out.aegis_proteins_main, aegis_short_reads.out.aegis_proteins_all)
      }
    } else {
      println "Skipping the Diamond2GO process because EDTA was not run (params.EDTA = '${params.EDTA}'). To enable EDTA, and consequently AEGIS and Diamond2GO, set EDTA = 'yes' in the nextflow.config file."
    }

  // Outputs
  emit:
    aegis_gff            = params.use_long_reads ? aegis_long_reads.out.aegis_gff : aegis_short_reads.out.aegis_gff
    aegis_pkl            = params.use_long_reads ? aegis_long_reads.out.aegis_pkl : aegis_short_reads.out.aegis_pkl
    aegis_proteins_main  = params.use_long_reads ? aegis_long_reads.out.aegis_proteins_main : aegis_short_reads.out.aegis_proteins_main
    aegis_proteins_all   = params.use_long_reads ? aegis_long_reads.out.aegis_proteins_all : aegis_short_reads.out.aegis_proteins_all
    diamond2go_output    = diamond2go.out
}
