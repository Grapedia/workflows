nextflow.enable.dsl=2

include { aegis_merge } from "../modules/aegis_merge"
include { aegis_liftoff_gene_ids } from "../modules/aegis_liftoff_ids"
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
    // Old-assembly annotation transferred onto the new assembly coordinates
    // (old gene IDs, new coordinates), and the true pre-transfer annotation
    // on the old assembly (for synteny conservation checks). Used to carry
    // old gene IDs over onto their corresponding new-assembly gene.
    liftoff_gff3
    previous_annotations

  main:

    def aegis_merge_script = file("${projectDir}/scripts/run_aegis_merge.sh")
    def apply_liftoff_ids_script = file("${projectDir}/scripts/apply_liftoff_gene_ids.py")

    vitvi_annotation = aegis_merge(
      masked_genome,
      mikado_gff3,
      aegis_merge_script
    )

    merged_annotation = aegis_liftoff_gene_ids(
      vitvi_annotation.aegis_gff,
      vitvi_annotation.aegis_proteins_all,
      vitvi_annotation.aegis_proteins_main,
      liftoff_gff3,
      previous_annotations,
      apply_liftoff_ids_script
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
    liftoff_gene_id_correspondence = merged_annotation.correspondence_tsv
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
