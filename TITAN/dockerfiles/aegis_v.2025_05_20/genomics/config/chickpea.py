#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from config.paths import root, GO_file

species="chickpea"
# assembly/annotation folder if any
assembly=""
species_path=f"{root}{species}/"
path = f"{species_path}{assembly}"

genome_files = {}
genome_files["chickpea"] = f"{path}chickpea_assembly/GCF_000331145.1_ASM33114v1_genomic.fasta"

annotation_files = {}
annotation_files["chickpea_on_chickpea"] = f"{path}chickpea_annotation/GCF_000331145.1_ASM33114v1_genomic.gff3"

annotation_transfer_files = {}

genome_dictionary = {}

for tag in genome_files:
    genome_dictionary[tag] = {}

# note that feature should match genome id property as understood by biopython
f_in = open(genome_files["chickpea"])
for line in f_in:
    if line.startswith(">"):
        feature = line[1:-1].split(" ")[0]
        if "chromosome" in line:
            # careful as this does not work genomes with more than 9 chromosomes
            number = "{:02d}".format(int(line.split("chromosome Ca")[1][0]))
            genome_dictionary["chickpea"][feature] = f"chr{number}"
f_in.close()
