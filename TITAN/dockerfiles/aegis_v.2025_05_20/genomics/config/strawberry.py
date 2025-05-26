#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from config.paths import root, GO_file

species="strawberry"
# assembly/annotation folder if any
assembly=""
species_path=f"{root}{species}/"
path = f"{species_path}{assembly}"

genome_files = {}
genome_files["ananassa"] = f"{path}fragaria_ananassa_genome/F_ana_Camarosa_6-28-17.fasta"
genome_files["vesca"] = f"{path}fragaria_vesca_genome/Fragaria_vesca_v4.0.a1.fasta"

# notice _on_ nomenclature
annotation_files = {}
annotation_files["ananassa_on_ananassa"] = f"{path}fragaria_ananassa_annotation/Fragaria_ananassa_v1.0.a2.genes_mod_tidy.gff3"
annotation_files["vesca_on_vesca"] = f"{path}fragaria_vesca_annotation/Fragaria_vesca_v4.0.a2.genes_tidy.gff3"

# notice _from_ and _on_ nomenclature
# leave empty if no liftoff or equivalent files
annotation_transfer_files = {}

functional_annotation = {}

genome_dictionary = {}

for tag in genome_files:
    genome_dictionary[tag] = {}