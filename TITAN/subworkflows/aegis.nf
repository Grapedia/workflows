nextflow.enable.dsl=2

include { aegis_merge } from "../modules/aegis_merge"
include { aegis_liftoff_gene_ids } from "../modules/aegis_liftoff_ids"
include { diamond2go } from "../modules/diamond2go"
include { eggnog_mapper } from "../modules/eggnog_mapper"
include { interproscan } from "../modules/interproscan"
include { mender_split_genes } from "../modules/mender_split_genes"

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
    // Inputs of the optional Mender split-gene merge (params.run_mender): unmasked assembly (to
    // translate merged CDS), protein evidence FASTAs (reference proteome) and long-read BAMs
    // (assets/empty_long_reads.bam without long reads).
    new_assembly
    protein_evidence_fastas
    long_read_bams

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

    // Mender merges split gene models before the functional annotation, so every downstream step
    // (function, primary_noTE, QC) sees the merged genes.
    if (params.run_mender) {
      mended = mender_split_genes(
        merged_annotation.aegis_gff,
        merged_annotation.aegis_proteins_all,
        merged_annotation.aegis_proteins_main,
        new_assembly,
        protein_evidence_fastas,
        long_read_bams,
        params.aegis_gene_id_prefix,
        [
          evalue: params.mender_evalue,
          exclude_os: params.mender_reference_exclude_os,
          flags: params.mender_flags,
          skip_flags: params.mender_skip_flags,
          include_review: params.mender_include_review
        ],
        file("${projectDir}/scripts/run_mender_split_genes.sh")
      )
      final_gff = mended.gff
      final_proteins_all = mended.proteins_all
      final_proteins_main = mended.proteins_main
    } else {
      final_gff = merged_annotation.aegis_gff
      final_proteins_all = merged_annotation.aegis_proteins_all
      final_proteins_main = merged_annotation.aegis_proteins_main
    }

    functional_annotation = diamond2go(
      final_proteins_all,
      final_proteins_main
    )

    eggnog_annotation = eggnog_mapper(
      final_proteins_all,
      final_proteins_main
    )

    interproscan_annotation = interproscan(
      final_proteins_all,
      final_proteins_main
    )

  emit:
    aegis_gff            = final_gff
    aegis_proteins_all   = final_proteins_all
    aegis_proteins_main  = final_proteins_main
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
