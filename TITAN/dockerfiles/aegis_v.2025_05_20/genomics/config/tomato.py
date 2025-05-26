#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from config.paths import root, GO_file

species="tomato"
# assembly/annotation folder if any
assembly=""
species_path=f"{root}{species}/"
path = f"{species_path}{assembly}"

genome_files = {}
genome_files["ITAG4"] = f"{path}Heinz_ITAG4_genome/ITAG4.fasta"

# notice _on_ nomenclature
annotation_files = {}
annotation_files["4.0_4.2_merge_on_ITAG4"] = f"{path}aegis_output/gffs/4.0_4.2_merge_on_ITAG4_aegis_fcounts_dapfit.gff3"

# notice _from_ and _on_ nomenclature
# leave empty if no liftoff or equivalent files
annotation_transfer_files = {}

functional_annotation = {}
functional_annotation["4.0_4.2_merge_on_ITAG4_mapman"] = f"{path}/Heinz_ITAG4_annotation/functional_annotation/MapMan.Mercator.tomato.gmt"
functional_annotation["4.0_4.2_merge_on_ITAG4_interpro"] = f"{path}/Heinz_ITAG4_annotation/functional_annotation/4.0_4.2_merge_aegis_on_ITAG4_all_proteins.interpro"
functional_annotation["4.0_4.2_merge_on_ITAG4_eggnog"] = f"{path}/Heinz_ITAG4_annotation/functional_annotation/out.emapper.annotations"

genome_dictionary = {}

for tag in genome_files:
    genome_dictionary[tag] = {}