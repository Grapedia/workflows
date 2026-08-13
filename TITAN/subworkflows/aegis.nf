nextflow.enable.dsl=2

include { aegis_merge } from "../modules/aegis_merge"
include { diamond2go } from "../modules/diamond2go"
include { eggnog_mapper } from "../modules/eggnog_mapper"
include { interproscan } from "../modules/interproscan"

workflow aegis {

  take:
    masked_genome
    // Already-consolidated, per-locus-scored annotation from mikado_pick,
    // built from every evidence source (liftoff, egapx, BRAKER, helixer,
    // and all transcript assemblies). AEGIS only renames (Vitvi prefix)
    // and tidies it here.
    mikado_gff3

  main:

    def aegis_merge_script = file("${projectDir}/scripts/run_aegis_merge.sh")

    merged_annotation = aegis_merge(
      masked_genome,
      mikado_gff3,
      aegis_merge_script
    )

    functional_annotation = diamond2go(
      merged_annotation.aegis_proteins_all,
      merged_annotation.aegis_proteins_main
    )

    eggnog_annotation = eggnog_mapper(
      merged_annotation.aegis_proteins_all,
      merged_annotation.aegis_proteins_main
    )

    interproscan_annotation = interproscan(
      merged_annotation.aegis_proteins_all,
      merged_annotation.aegis_proteins_main
    )

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
    interproscan_all_tsv    = interproscan_annotation.proteins_all_tsv
    interproscan_main_tsv   = interproscan_annotation.proteins_main_tsv
    interproscan_versions   = interproscan_annotation.versions
}
