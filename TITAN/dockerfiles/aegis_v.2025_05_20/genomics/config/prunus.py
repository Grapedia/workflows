#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from config.paths import root, GO_file

species="prunus"
# assembly/annotation folder if any
assembly=""
species_path=f"{root}{species}/"
path = f"{species_path}{assembly}"

genome_files = {}
genome_files["persica"] = f"{path}persica_v2_assembly/Prunus_persica_v2.0.a1_scaffolds.fasta"
genome_files["davidiana"] = f"{path}davidiana_assembly/Pda.genome.update.fa"
genome_files["dulcis_f0"] = f"{path}dulcis_assembly/Texas_F0_K80_chr.fasta"
genome_files["dulcis_f1"] = f"{path}dulcis_assembly/Texas_F1_K80_chr.fasta"

# notice _on_ nomenclature
annotation_files = {}
annotation_files["persica_on_persica"] = f"{path}persica_v2_annotation/Prunus_persica_v2.0.a1.gene.gff3"
annotation_files["davidiana_on_davidiana"] = f"{path}davidiana_annotation/Pda.update.gff3"
annotation_files["dulcis_f0_on_dulcis_f0"] = f"{path}dulcis_annotation/Texas_F0_gene_models.gff3"
annotation_files["dulcis_f1_on_dulcis_f1"] = f"{path}dulcis_annotation/Texas_F1_gene_models.gff3"

# notice _from_ and _on_ nomenclature
# leave empty if no liftoff or equivalent files
annotation_transfer_files = {}
annotation_files["davidiana_flanking_on_persica"] = f"{path}persica_v2_annotation/transferred_annotations/Pdavidiana_liftoff_to_Ppersica_copies_flanking_0.2.gff3"
annotation_files["davidiana_on_persica"] = f"{path}persica_v2_annotation/transferred_annotations/Pdavidiana_liftoff_to_Ppersica_copies.gff3"
annotation_files["dulcis_f1_on_persica"] = f"{path}persica_v2_annotation/transferred_annotations/Pdulcis_liftoff_to_Ppersica_no_copies.gff3"
annotation_files["persica_on_dulcis_f1"] = f"{path}dulcis_annotation/transferred_annotations/Ppersica_liftoff_to_Pdulcis_no_copies.gff3"
annotation_files["persica_on_davidiana"] = f"{path}davidiana_annotation/transferred_annotations/Ppersica_liftoff_to_Pdavidiana_copies.gff3"

functional_annotation = {}

genome_dictionary = {}

for tag in genome_files:
    genome_dictionary[tag] = {}

f_in = open(genome_files["davidiana"])
for line in f_in:
    if line.startswith(">"):
        feature = line[1:].strip()
        if feature.startswith("Pda"):
            number = "{:02d}".format(int(line.split("Pda")[-1]))
            genome_dictionary["davidiana"][feature] = f"chr{number}"
f_in.close()

f_in = open(genome_files["persica"])
for line in f_in:
    if line.startswith(">"):
        feature = line[1:].strip()
        if feature.startswith("Pp"):
            number = "{:02d}".format(int(line.split("Pp")[-1]))
            genome_dictionary["persica"][feature] = f"chr{number}"
f_in.close()