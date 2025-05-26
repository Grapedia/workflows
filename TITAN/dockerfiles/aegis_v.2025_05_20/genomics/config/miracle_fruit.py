#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from config.paths import root, GO_file

species="miracle_fruit"
# assembly/annotation folder if any
assembly=""
species_path=f"{root}{species}/"
path = f"{species_path}{assembly}"

genome_files = {}
genome_files["miracle_fruit"] = f"{path}miracle_fruit_genome/miracle_fruit.fasta"

# notice _on_ nomenclature
annotation_files = {}
annotation_files["miracle_fruit_on_miracle_fruit"] = f"{path}miracle_fruit_annotation/miracle_fruit_tidy.gff3"

# notice _from_ and _on_ nomenclature
# leave empty if no liftoff or equivalent files
annotation_transfer_files = {}

functional_annotation = {}

genome_dictionary = {}

for tag in genome_files:
    genome_dictionary[tag] = {}
