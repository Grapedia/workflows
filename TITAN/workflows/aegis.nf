nextflow.enable.dsl=2

// Include various processing modules
include { aegis_short_reads } from "../modules/aegis_short_reads"
include { aegis_long_reads } from "../modules/aegis_long_reads"
include { diamond2go } from "../modules/diamond2go"

workflow aegis {

  take:
    input_data

  main:

    if (params.EDTA == 'yes') {
      masked_genome = input_data.out.get("masked_genome")
      println "Masked genome file: ${masked_genome}"
    }

    if (params.use_long_reads) {
      Stringtie_merging_long_reads_default_args = input_data.out.get("Stringtie_merging_long_reads_default_args")
      println "Stringtie_merging_long_reads_default_args file: ${Stringtie_merging_long_reads_default_args}"

      Stringtie_merging_long_reads_alt_args = input_data.out.get("Stringtie_merging_long_reads_alt_args")
      println "Stringtie_merging_long_reads_alt_args file: ${Stringtie_merging_long_reads_alt_args}"
    }
    
    braker3_prediction_augustus = input_data.out.get("braker3_augustus")
    println "braker3_prediction_augustus file: ${braker3_prediction_augustus}"

    braker3_prediction_genemark = input_data.out.get("braker3_genemark")
    println "braker3_prediction_genemark file: ${braker3_prediction_genemark}"

    liftoff_annotation = input_data.out.get("liftoff_annotations")
    println "liftoff_annotation file: ${liftoff_annotation}"

    Stringtie_merging_short_reads_STAR_default_args_stranded = input_data.out.get("stringtie_short_reads_STAR_default_stranded")
    println "Stringtie_merging_short_reads_STAR_default_args_stranded file: ${Stringtie_merging_short_reads_STAR_default_args_stranded}"

    Stringtie_merging_short_reads_STAR_alt_args_stranded = input_data.out.get("stringtie_short_reads_STAR_alt_stranded")
    println "Stringtie_merging_short_reads_STAR_alt_args_stranded file: ${Stringtie_merging_short_reads_STAR_alt_args_stranded}"

    gffcompare_star_psiclass_stranded = input_data.out.get("gffcompare_star_psiclass_stranded")
    println "gffcompare_star_psiclass_stranded file: ${gffcompare_star_psiclass_stranded}"

    Stringtie_merging_short_reads_STAR_default_args_unstranded = input_data.out.get("stringtie_short_reads_STAR_default_unstranded")
    println "Stringtie_merging_short_reads_STAR_default_args_unstranded file: ${Stringtie_merging_short_reads_STAR_default_args_unstranded}"

    Stringtie_merging_short_reads_STAR_alt_args_unstranded = input_data.out.get("stringtie_short_reads_STAR_alt_unstranded")
    println "Stringtie_merging_short_reads_STAR_alt_args_unstranded file: ${Stringtie_merging_short_reads_STAR_alt_args_unstranded}"

    gffcompare_star_psiclass_unstranded = input_data.out.get("gffcompare_star_psiclass_unstranded")
    println "gffcompare_star_psiclass_unstranded file: ${gffcompare_star_psiclass_unstranded}"

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
        masked_genome,
        braker3_prediction_with_long_reads_augustus,
        braker3_prediction_with_long_reads_genemark,
        liftoff_annotation,
        Stringtie_merging_long_reads_default_args,
        Stringtie_merging_long_reads_alt_args,
        Stringtie_merging_short_reads_STAR_default_args_stranded,
        Stringtie_merging_short_reads_STAR_alt_args_stranded,
        gffcompare_star_psiclass_stranded,
        gffcompare_star_psiclass_unstranded,
        Stringtie_merging_short_reads_STAR_default_args_unstranded,
        Stringtie_merging_short_reads_STAR_alt_args_unstranded
        )
      } else {
        aegis_short_reads(
        file(params.new_assembly).getParent(),
        file(params.new_assembly).getName(),
        file(params.protein_samplesheet).getParent(),
        file(params.protein_samplesheet).getName(),
        masked_genome,
        braker3_prediction_augustus,
        braker3_prediction_genemark,
        liftoff_annotation,
        Stringtie_merging_short_reads_STAR_default_args_stranded,
        Stringtie_merging_short_reads_STAR_alt_args_stranded,
        gffcompare_star_psiclass_stranded,
        gffcompare_star_psiclass_unstranded,
        Stringtie_merging_short_reads_STAR_default_args_unstranded,
        Stringtie_merging_short_reads_STAR_alt_args_unstranded
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

    // Define a dictionnary (Map) to store outputs for aegis workflow
    def outputs_map = [:]

    if (params.use_long_reads) {
      outputs_map["aegis_gff"] = aegis_long_reads.out.aegis_gff.publishDir(output_directory)
      outputs_map["aegis_pkl"] = aegis_long_reads.out.aegis_pkl.publishDir(output_directory)
      outputs_map["aegis_proteins_main"] = aegis_long_reads.out.aegis_proteins_main.publishDir(output_directory)
      outputs_map["aegis_proteins_all"] = aegis_long_reads.out.aegis_proteins_all.publishDir(output_directory)
    } else {
      outputs_map["aegis_gff"] = aegis_short_reads.out.aegis_gff.publishDir(output_directory)
      outputs_map["aegis_pkl"] = aegis_short_reads.out.aegis_pkl.publishDir(output_directory)
      outputs_map["aegis_proteins_main"] = aegis_short_reads.out.aegis_proteins_main.publishDir(output_directory)
      outputs_map["aegis_proteins_all"] = aegis_short_reads.out.aegis_proteins_all.publishDir(output_directory)
    }

    // Results of gffcompare after PsiClass for stranded and unstranded data
    outputs_map["diamond2go"] = diamond2go.out.publishDir(output_directory)

  // Outputs
  emit: outputs_map
}
