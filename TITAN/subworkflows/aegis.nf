nextflow.enable.dsl=2

// Include various processing modules
include { aegis_merge } from "../modules/aegis_merge"
include { diamond2go } from "../modules/diamond2go"
include { eggnog_mapper } from "../modules/eggnog_mapper"

workflow aegis {

  take:
    // Named evidence files consumed by AEGIS. Long-read GTF inputs may be empty
    // sentinel files when no long-read branch is active; the short-read module
    // ignores them by construction.
    masked_genome
    braker_augustus_gff
    braker_genemark_gtf
    liftoff_annotation
    egapx_gff3
    star_stringtie_default_stranded
    star_stringtie_alt_stranded
    star_psiclass_stranded
    star_psiclass_unstranded
    star_stringtie_default_unstranded
    star_stringtie_alt_unstranded
    long_reads_default
    long_reads_alt
    has_long_reads

  main:

    def aegis_mode = has_long_reads ? 'short_and_long_reads' : 'short_reads'
    def aegis_merge_script = file("${projectDir}/scripts/run_aegis_merge.sh")

    // ----------------------------------------------------------------------------------------
    //     AEGIS CLI merge over all named annotation evidence
    // ----------------------------------------------------------------------------------------

    merged_annotation = aegis_merge(
      aegis_mode,
      masked_genome,
      braker_augustus_gff,
      braker_genemark_gtf,
      liftoff_annotation,
      egapx_gff3,
      long_reads_default,
      long_reads_alt,
      star_stringtie_default_stranded,
      star_stringtie_alt_stranded,
      star_psiclass_stranded,
      star_psiclass_unstranded,
      star_stringtie_default_unstranded,
      star_stringtie_alt_unstranded,
      aegis_merge_script
    )

    // ----------------------------------------------------------------------------------------
    //               Diamond2GO on proteins predicted with TITAN/Aegis
    // ----------------------------------------------------------------------------------------

    functional_annotation = diamond2go(
      merged_annotation.aegis_proteins_all,
      merged_annotation.aegis_proteins_main
    )

    // ----------------------------------------------------------------------------------------
    //     eggNOG-mapper on the same AEGIS-derived proteins
    // ----------------------------------------------------------------------------------------
    // Always invoked, like diamond2go; the module itself emits placeholder
    // outputs and skips emapper.py when params.run_eggnog_mapper is false.

    eggnog_annotation = eggnog_mapper(
      merged_annotation.aegis_proteins_all,
      merged_annotation.aegis_proteins_main
    )

  // Outputs
  emit:
    aegis_gff            = merged_annotation.aegis_gff
    aegis_proteins_all   = merged_annotation.aegis_proteins_all
    aegis_proteins_main  = merged_annotation.aegis_proteins_main
    aegis_versions       = merged_annotation.versions
    diamond2go_all       = functional_annotation.proteins_all_diamond
    diamond2go_main      = functional_annotation.proteins_main_diamond
    diamond2go_versions  = functional_annotation.versions
    eggnog_annotations_all  = eggnog_annotation.proteins_all_annotations
    eggnog_annotations_main = eggnog_annotation.proteins_main_annotations
    eggnog_versions         = eggnog_annotation.versions
}
