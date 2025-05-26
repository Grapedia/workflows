#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from config.paths import root, GO_file

species="brachypodium"
# assembly/annotation folder if any
assembly=""
species_path=f"{root}{species}/"
path = f"{species_path}{assembly}"

genome_files = {}
genome_files["bdistachyon"] = f"{path}bdistachyon_genome/Bdistachyon_556_v3.0.fa"

# notice _on_ nomenclature
annotation_files = {}
annotation_files["bdistachyon_on_bdistachyon"] = f"{path}bdistachyon_annotation/Bdistachyon_556_v3.2.gene_exons.gff3"

# notice _from_ and _on_ nomenclature
# leave empty if no liftoff or equivalent files
annotation_transfer_files = {}

functional_annotation = {}

genome_dictionary = {}

for tag in genome_files:
    genome_dictionary[tag] = {}

# note that feature should match genome id property as understood by biopython
f_in = open(genome_files["bdistachyon"])
for line in f_in:
    if line.startswith(">"):
        feature = line[1:-1].split(" ")[0]
        if line.startswith(">Bd"):
            number_string = line.split("Bd")[1]
            if len(number_string) < 4:
                number = "{:02d}".format(int(number_string))
                genome_dictionary["bdistachyon"][feature] = f"chr{number}"
f_in.close()