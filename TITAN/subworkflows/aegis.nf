nextflow.enable.dsl=2

// Include various processing modules
include { aegis_short_reads } from "../modules/aegis_short_reads"
include { aegis_long_reads } from "../modules/aegis_long_reads"
include { diamond2go } from "../modules/diamond2go"

workflow aegis {

  take:
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

    def has_long_reads_enabled = has_long_reads

    // ----------------------------------------------------------------------------------------
    //     AEGIS CLI merge over all named annotation evidence
    // ----------------------------------------------------------------------------------------

    if (has_long_reads_enabled) {
      aegis_long_reads(
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
        star_stringtie_alt_unstranded
        )
    } else {
      aegis_short_reads(
        masked_genome,
        braker_augustus_gff,
        braker_genemark_gtf,
        liftoff_annotation,
        egapx_gff3,
        star_stringtie_default_stranded,
        star_stringtie_alt_stranded,
        star_psiclass_stranded,
        star_psiclass_unstranded,
        star_stringtie_default_unstranded,
        star_stringtie_alt_unstranded
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
