#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from config.paths import root, GO_file

species="algae"
# assembly/annotation folder if any
assembly="chlorella_vulgaris/"
species_path=f"{root}{species}/"
path = f"{species_path}{assembly}"

genome_files = {}
genome_files["cvul"] = f"{path}cvul_assembly/GCA_023343905.1_cvul_genomic.fasta"

# notice _on_ nomenclature
annotation_files = {}
annotation_files["cvul_on_cvul"] = f"{path}cvul_annotation/genomic.gff"

# notice _from_ and _on_ nomenclature
# leave empty if no liftoff or equivalent files
annotation_transfer_files = {}

functional_annotation = {}

genome_dictionary = {}

for tag in genome_files:
    genome_dictionary[tag] = {}
