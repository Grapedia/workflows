#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from config.paths import root, GO_file

species="hop"
# assembly/annotation folder if any
assembly=""
species_path=f"{root}{species}/"
path = f"{species_path}{assembly}"

genome_files = {}
genome_files["dovetail_cascade"] = f"{path}dovetail_cascade_genome/hop_unmasked_genome.fasta"

# notice _on_ nomenclature
annotation_files = {}
annotation_files["hop_merge_on_dovetail_cascade"] = f"{path}dovetail_cascade_annotation/hop_annotation_tidy.gff3"

# notice _from_ and _on_ nomenclature
# leave empty if no liftoff or equivalent files
annotation_transfer_files = {}

functional_annotation = {}

genome_dictionary = {}

for tag in genome_files:
    genome_dictionary[tag] = {}