#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from config.paths import root, GO_file

species="algae"
# assembly/annotation folder if any
assembly="chlamydomonas_reinhardtii/"
species_path=f"{root}{species}/"
path = f"{species_path}{assembly}"

genome_files = {}
genome_files["Chlre_5.6"] = f"{path}5.6_assembly/Chlre5_6_AssemblyScaffolds.fasta"
genome_files["Chlre_6.0"] = f"{path}6.0_assembly/CreinhardtiiCC_4532_707_v6.0.fasta"

# notice _on_ nomenclature
annotation_files = {}
annotation_files["Chlre_5.6_on_Chlre_5.6"] = f"{path}5.6_annotation/Chlre5_6_GeneCatalog_20200117_tidy.gff3"
annotation_files["Chlre_6.1_on_Chlre_6.0"] = f"{path}6.1_annotation/CreinhardtiiCC_4532_707_v6.1.gene_exons_tidy.gff3"


functional_annotation = {}
functional_annotation["Chlre_6.1_on_Chlre_6.0_mapman"] = f"{path}/6.1_annotation/functional_annotation/Chlre_6.1_on_Chlre_6.0_mapman.gmt"
functional_annotation["Chlre_6.1_on_Chlre_6.0_interpro"] = f"{path}/6.1_annotation/functional_annotation/Chlre_6.1_on_Chlre_6.0_all_proteins.interpro"
functional_annotation["Chlre_6.1_on_Chlre_6.0_eggnog"] = f"{path}/6.1_annotation/functional_annotation/out.emapper.annotations"

# notice _from_ and _on_ nomenclature
# leave empty if no liftoff or equivalent files
annotation_transfer_files = {}

genome_dictionary = {}

for tag in genome_files:
    genome_dictionary[tag] = {}

