#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from config.paths import root, GO_file

species="grapevine"
# assembly/annotation folder if any
assembly="PN40024/"
species_path=f"{root}{species}/"
path = f"{species_path}{assembly}"

genome_files = {}
genome_files["8x"] = f"{path}8X_genome_genoscope/8x.fasta"
genome_files["12XNCBI"] = f"{path}12X_genome_NCBI/12X_genome_NCBI.fasta"
genome_files["12Xv2"] = f"{path}12Xv2_genome/12Xv2_genome.fasta"
genome_files["40X_ref"] = f"{path}40X_genome/v4_genome_ref.fasta"
genome_files["40X_alt"] = f"{path}40X_genome/v4_genome_alt.fasta"
genome_files["40X_ref_hard"] = f"{path}40X_genome/v4_genome_ref_hard_mod1.fasta"
genome_files["T2T_ref"] = f"{path}T2T_genome/T2T_ref.fasta"
genome_files["T2T_alt"] = f"{path}T2T_genome/T2T_alt.fasta"


pre_annotation_files = {}
pre_annotation_files["just_lncRNAs_on_T2T_ref"] = f"{path}T2T_annotation/POTENTIAL_LNCRNAS_pred.gff3"
pre_annotation_files["v5.1_ncRNAs_removed_on_T2T_ref"] = f"{path}T2T_annotation/v5.1_ncRNAs_removed.gff3"

# notice _on_ nomenclature
annotation_files = {}
annotation_files["8x_on_8x"] = f"{path}8X_annotation_genoscope/8x_tidy.gff3"
annotation_files["12XNCBI_on_12XNCBI"] = f"{path}12X_annotation_NCBI/NCBI_original_tidy_loc.gff3"
annotation_files["v0_on_12Xv2"] = f"{path}12Xv2_annotation_transfer_URGI/v0_on_12Xv2_tidy.gff3"
annotation_files["v1_on_12Xv2"] = f"{path}12Xv2_annotation_transfer_URGI/v1_on_12Xv2_tidy.gff3"
annotation_files["v2_on_12Xv2"] = f"{path}12Xv2_annotation_transfer_URGI/v2_on_12Xv2_tidy.gff3"
annotation_files["v3_on_12Xv2"] = f"{path}12Xv2_annotation/v3_tidy.gff3"
annotation_files["cyps_on_12Xv2"] = f"{path}12Xv2_annotation/cyps_camille_paper_tidy_fc_simple.gff3"
annotation_files["cyps_renamed_on_12Xv2"] = f"{path}12Xv2_annotation/cyps_renamed_on_12Xv2_aegis_dapfit.gff3"
annotation_files["v4.1_on_40X_ref"] = f"{path}40X_annotation/v4_1_ref_tidy.gff3"
annotation_files["v4.1_on_40X_alt"] = f"{path}40X_annotation/v4_1_alt_tidy.gff3"
annotation_files["v4.3_on_40X_ref"] = f"{path}40X_annotation/v4_3_ref_tidy.gff3"
annotation_files["5.0_on_T2T_ref"] = f"{path}T2T_annotation/5.0_on_T2T_ref.gff3"
annotation_files["5.1_on_T2T_ref"] = f"{path}T2T_annotation/5.1_on_T2T_ref.gff3"

#test_files
# annotation_files["transcriptomes_on_40X_ref"] = f"{path}40X_annotation_tests/merged_transcriptomes_without_chloroplast_mitochondrion.gff3"
# annotation_files["geneid_on_40X_ref"] = f"{path}40X_annotation_tests/geneid.gff"
# annotation_files["glimmer_on_40X_ref"] = f"{path}40X_annotation_tests/glimmerhmm.gff"
# annotation_files["augustus_on_40X_ref"] = f"{path}40X_annotation_tests/augustus.gff"
# annotation_files["v3_from_12Xv2_raw_on_40X_ref"] = f"{path}40X_annotation_transfer/original/v3_copies_liftoff_ref_raw.gff3"
# annotation_files["exonerate_on_40X_ref"] = f"{path}40X_annotation_tests/merged_protein_final_alignments.gff3"
# annotation_files["arabidopsis_on_40X_ref"] = f"{path}40X_annotation_tests/Arabidopsis.gff3"

# annotation_files["genemark_on_40X_ref"] = f"{path}40X_annotation_tests/genemark.gff"
# annotation_files["final_reduced_on_40X_ref"] = f"{path}40X_annotation_tests/final_annotation_test36_5.gff3"
# annotation_files["final_reduced_without_PN40024_on_40X_ref"] = f"{path}40X_annotation_tests/final_annotation_test36_noPN.gff3"
# #annotation_files["EVM_on_40X_ref"] = f"{path}40X_annotation_tests/EVM_unfiltered_V4_test.gff3"
# annotation_files["EVM_PASA_on_40X_ref"] = f"{path}40X_annotation_tests/pasa.gene_structures_post_PASA_updates.1948239.gff3"

# annotation_files["psiclass_stranded_on_40X_ref"] = f"{path}40X_annotation_tests/psiclass_stranded_STAR.gff"
# annotation_files["psiclass_unstranded_on_40X_ref"] = f"{path}40X_annotation_tests/psiclass_unstranded_STAR.gff"
# annotation_files["stringtie_default_on_40X_ref"] = f"{path}40X_annotation_tests/stringtie_SRA1_default_STAR_2.gff"
# annotation_files["stringtie_morus_on_40X_ref"] = f"{path}40X_annotation_tests/stringtie_SRA1_MorusCommands_STAR_2.gff"
# annotation_files["stringtie_stranded_default_on_40X_ref"] = f"{path}40X_annotation_tests/stringtie_stranded_default_STAR_2.gff"
# annotation_files["stringtie_stranded_morus_on_40X_ref"] = f"{path}40X_annotation_tests/stringtie_stranded_MorusCommands_STAR_2.gff"
# annotation_files["stringtie_unstranded_on_40X_ref"] = f"{path}40X_annotation_tests/stringtie_unstranded_default_STAR_2.gff"
# annotation_files["stringtie_unstranded_morus_on_40X_ref"] = f"{path}40X_annotation_tests/stringtie_unstranded_MorusCommands_STAR_2.gff"

# notice _from_ and _on_ nomenclature
# leave empty if no liftoff or equivalent files
annotation_transfer_files = {}
annotation_transfer_files["8x_from_8x_on_40X_ref"] = f"{path}40X_annotation_transfer/8x_copies_liftoff_ref_mod.gff3"
annotation_transfer_files["8x_from_8x_on_40X_alt"] = f"{path}40X_annotation_transfer/8x_copies_liftoff_alt_mod.gff3"
#annotation_transfer_files["12XNCBI_from_12XNCBI_on_40X_ref"] = f"{path}40X_annotation_transfer/NCBI_original_copies_liftoff_ref.gff3"
#annotation_transfer_files["12XNCBI_from_12XNCBI_on_40X_alt"] = f"{path}40X_annotation_transfer/NCBI_original_copies_liftoff_alt.gff3"
annotation_transfer_files["v0_from_12Xv2_on_40X_ref"] = f"{path}40X_annotation_transfer/v0_on_12Xv2_copies_liftoff_ref.gff3"
annotation_transfer_files["v0_from_12Xv2_on_40X_alt"] = f"{path}40X_annotation_transfer/v0_on_12Xv2_copies_liftoff_alt.gff3"
annotation_transfer_files["v1_from_12Xv2_on_40X_ref"] = f"{path}40X_annotation_transfer/v1_on_12Xv2_copies_liftoff_ref.gff3"
annotation_transfer_files["v1_from_12Xv2_on_40X_alt"] = f"{path}40X_annotation_transfer/v1_on_12Xv2_copies_liftoff_alt.gff3"
annotation_transfer_files["v2_from_12Xv2_on_40X_ref"] = f"{path}40X_annotation_transfer/v2_on_12Xv2_copies_liftoff_ref.gff3"
annotation_transfer_files["v2_from_12Xv2_on_40X_alt"] = f"{path}40X_annotation_transfer/v2_on_12Xv2_copies_liftoff_alt.gff3"
annotation_transfer_files["v3_from_12Xv2_on_40X_ref"] = f"{path}40X_annotation_transfer/v3_copies_liftoff_ref_clean.gff3"
annotation_transfer_files["v3_from_12Xv2_on_40X_alt"] = f"{path}40X_annotation_transfer/v3_copies_liftoff_alt_clean.gff3"

annotation_transfer_files["8x_from_8x_on_T2T_ref"] = f"{path}T2T_annotation_transfer/8x_copies_liftoff_T2T_ref.gff3"
#annotation_transfer_files["8x_from_8x_on_T2T_alt"] = f"{path}T2T_annotation_transfer/8x_copies_liftoff_T2T_alt.gff3"
annotation_transfer_files["12XNCBI_from_12XNCBI_on_T2T_ref"] = f"{path}T2T_annotation_transfer/NCBI_original_copies_liftoff_T2T_ref.gff3"
#annotation_transfer_files["12XNCBI_from_12XNCBI_on_T2T_alt"] = f"{path}T2T_annotation_transfer/NCBI_original_copies_liftoff_T2T_alt.gff3"
annotation_transfer_files["v0_from_12Xv2_on_T2T_ref"] = f"{path}T2T_annotation_transfer/v0_on_12Xv2_copies_liftoff_T2T_ref.gff3"
#annotation_transfer_files["v0_from_12Xv2_on_T2T_alt"] = f"{path}T2T_annotation_transfer/v0_on_12Xv2_copies_liftoff_T2T_alt.gff3"
annotation_transfer_files["v1_from_12Xv2_on_T2T_ref"] = f"{path}T2T_annotation_transfer/v1_on_12Xv2_copies_liftoff_T2T_ref.gff3"
#annotation_transfer_files["v1_from_12Xv2_on_T2T_alt"] = f"{path}T2T_annotation_transfer/v1_on_12Xv2_copies_liftoff_T2T_alt.gff3"
annotation_transfer_files["v2_from_12Xv2_on_T2T_ref"] = f"{path}T2T_annotation_transfer/v2_on_12Xv2_copies_liftoff_T2T_ref.gff3"
#annotation_transfer_files["v2_from_12Xv2_on_T2T_alt"] = f"{path}T2T_annotation_transfer/v2_on_12Xv2_copies_liftoff_T2T_alt.gff3"
annotation_transfer_files["v3_from_12Xv2_on_T2T_ref"] = f"{path}T2T_annotation_transfer/v3_copies_liftoff_T2T_ref_clean.gff3"
#annotation_transfer_files["v3_from_12Xv2_on_T2T_alt"] = f"{path}T2T_annotation_transfer/v3_copies_liftoff_T2T_alt_clean.gff3"
annotation_transfer_files["cyps_renamed_from_12Xv2_on_T2T_ref"] = f"{path}T2T_annotation_transfer/cyps_renamed_default_liftoff_T2T_ref.gff3"
#annotation_transfer_files["cyps_renamed_from_12Xv2_on_T2T_alt"] = f"{path}T2T_annotation_transfer/cyps_renamed_default_liftoff_T2T_alt.gff3"
annotation_transfer_files["v4.3_from_40X_ref_on_T2T_ref"] = f"{path}T2T_annotation_transfer/v4_3_ref_copies_liftoff_T2T_ref.gff3"

functional_annotation = {}
functional_annotation["5.1_on_T2T_ref_mapman"] = f"{path}/T2T_annotation/functional_annotation/MapMan.Mercator.vitisv5.1aegis.gmt"
functional_annotation["5.1_on_T2T_ref_interpro"] = f"{path}/T2T_annotation/functional_annotation/5.1_aegis_on_T2T_ref_all_proteins.interpro"
functional_annotation["5.1_on_T2T_ref_eggnog"] = f"{path}/T2T_annotation/functional_annotation/out.emapper.annotations"

functional_annotation["v3_on_12Xv2_mapman"] = f"{path}/12Xv2_annotation/functional_annotation/Ziva_onto_with_antho_corrected_without_glucosinolates_PME_PG_PL.gmt"
functional_annotation["v3_on_12Xv2_interpro"] = f"{path}/12Xv2_annotation/functional_annotation/v3_aegis_on_12Xv2_all_proteins.interpro"
functional_annotation["v3_on_12Xv2_eggnog"] = f"{path}/12Xv2_annotation/functional_annotation/out.emapper.annotations"

symbols = {}
symbols["5.1_on_T2T_ref_symbols"] = f"{path}/T2T_annotation/functional_annotation/MapMan.Mercator.vitisv5.1aegis.gmt"

genome_dictionary = {}

for tag in genome_files:
    genome_dictionary[tag] = {}