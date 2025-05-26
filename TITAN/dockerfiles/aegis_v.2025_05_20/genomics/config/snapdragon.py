#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from config.paths import root, GO_file

from config.thale_cress import genome_files

species="snapdragon"
# assembly/annotation folder if any
assembly="ji7_assembly_3/"
species_path=f"{root}{species}/"
path = f"{species_path}{assembly}"

#genome_files = {}
genome_files["ji7_3"] = f"{path}ji7_genome/ama_ji7_genome.fasta"

# notice _on_ nomenclature
annotation_files = {}
annotation_files["ji7_3_on_ji7_3"] = f"{path}ji7_annotation/ama_ji7_annotation_tidy.gff3"

# notice _from_ and _on_ nomenclature
# leave empty if no liftoff or equivalent files
annotation_transfer_files = {}
annotation_transfer_files["Araport11_from_TAIR10_on_ji7_3"] = f"{path}ji7_annotation_transfer/Araport11_copies_liftoff_ji7_3_raw.gff3"

functional_annotation = {}
functional_annotation["ji7_3_on_ji7_3_mapman"] = f"{path}/ji7_annotation/functional_annotation/MapMan.Mercator.Snapdragon.gmt"
functional_annotation["ji7_3_on_ji7_3_interpro"] = f"{path}/ji7_annotation/functional_annotation/ji7_3_on_ji7_3_main_proteins_interpro.csv"
functional_annotation["ji7_3_on_ji7_3_eggnog"] = f"{path}/ji7_annotation/functional_annotation/out.emapper.annotations"

genome_dictionary = {}

for tag in genome_files:
    genome_dictionary[tag] = {}