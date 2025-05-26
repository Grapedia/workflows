#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from config.paths import root, GO_file

species="pisum"
# assembly/annotation folder if any
assembly=""
species_path=f"{root}{species}/"
path = f"{species_path}{assembly}"

genome_files = {}
genome_files["psativum"] = f"{path}psativum_assembly/pisum_sativum.fa"

# notice _on_ nomenclature
annotation_files = {}
annotation_files["psativum_on_psativum"] = f"{path}psativum_annotation/pisum_sativum.gff3"

# notice _from_ and _on_ nomenclature
# leave empty if no liftoff or equivalent files
annotation_transfer_files = {}

functional_annotation = {}

genome_dictionary = {}

for tag in genome_files:
    genome_dictionary[tag] = {}

# note that feature should match genome id property as understood by biopython
f_in = open(genome_files["psativum"])
for line in f_in:
    if line.startswith(">"):
        feature = line[1:-1].split(" ")[0]
        if "chromosome" in line:
            # careful as this does not work genomes with more than 9 chromosomes
            number = "{:02d}".format(int(line.split("Zhongwan6 chromosome ")[1][0]))
            genome_dictionary["psativum"][feature] = f"chr{number}"
        if "mitochondrion" in line:
            genome_dictionary["psativum"][feature] = f"chrM"
        elif "chloroplast" in line:
            genome_dictionary["psativum"][feature] = f"chrC"

f_in.close()

