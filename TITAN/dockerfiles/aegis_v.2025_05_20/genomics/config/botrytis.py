#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from config.paths import root, GO_file

species="botrytis"
# assembly/annotation folder if any
assembly=""
species_path=f"{root}{species}/"
path = f"{species_path}{assembly}"

genome_files = {}
genome_files["botrytis"] = f"{path}botrytis_cinerea_assembly/GCA_000349525.1_BcDW1_genomic.fasta"

annotation_files = {}
annotation_files["botrytis_on_botrytis"] = f"{path}botrytis_cinerea_annotation/Bcinerea_BcDW1.GCA_000349525.1.gff3"

annotation_transfer_files = {}

genome_dictionary = {}

for tag in genome_files:
    genome_dictionary[tag] = {}
