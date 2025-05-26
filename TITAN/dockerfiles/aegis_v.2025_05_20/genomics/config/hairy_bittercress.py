#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from config.paths import root, GO_file

species="hairy_bittercress"
# assembly/annotation folder if any
assembly=""
species_path=f"{root}{species}/"
path = f"{species_path}{assembly}"

genome_files = {}
genome_files["chi_v1"] = f"{path}cardamine_hirsuta_genome/chi_v1.fasta"

# notice _on_ nomenclature
annotation_files = {}
annotation_files["carhr38_on_chi_v1"] = f"{path}cardamine_hirsuta_annotation/carhr38_tidy.gff3"

# notice _from_ and _on_ nomenclature
# leave empty if no liftoff or equivalent files
annotation_transfer_files = {}

functional_annotation = {}

genome_dictionary = {}

for tag in genome_files:
    genome_dictionary[tag] = {}
