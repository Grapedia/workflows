#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from config.paths import root, GO_file

species="mulberry"
# assembly/annotation folder if any
assembly=""
species_path=f"{root}{species}/"
path = f"{species_path}{assembly}"

genome_files = {}
genome_files["Mal_HE_154X"] = f"{path}Mal_HE_154X_genome/Mal_HE_154X_genome.fasta"

# notice _on_ nomenclature
annotation_files = {}
annotation_files["Mal_HE_154X_1_on_Mal_HE_154X"] = f"{path}Mal_HE_154X_annotation/Mal_HE_154X_1_tidy.gff3"

# notice _from_ and _on_ nomenclature
# leave empty if no liftoff or equivalent files
annotation_transfer_files = {}

functional_annotation = {}
functional_annotation["Mal_HE_154X_1_on_Mal_HE_154X_mapman"] = f"{path}/Mal_HE_154X_annotation/functional_annotation/MapMan.Mal_HE_154X_1.v2.gmt"
functional_annotation["Mal_HE_154X_1_on_Mal_HE_154X_interpro"] = f"{path}/Mal_HE_154X_annotation/functional_annotation/Mal_HE_154X_1_aegis_on_Mal_HE_154X_all_proteins.interpro"
functional_annotation["Mal_HE_154X_1_on_Mal_HE_154X_eggnog"] = f"{path}/Mal_HE_154X_annotation/functional_annotation/out.emapper.annotations"

genome_dictionary = {}

for tag in genome_files:
    genome_dictionary[tag] = {}