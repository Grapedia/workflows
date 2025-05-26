#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from config.paths import root, GO_file

species="thale_cress"
# assembly/annotation folder if any
assembly=""
species_path=f"{root}{species}/"
path = f"{species_path}{assembly}"

genome_files = {}
genome_files["TAIR10"] = f"{path}TAIR10_genome/TAIR10_genome.fasta"

# notice _on_ nomenclature
annotation_files = {}
annotation_files["Araport11_on_TAIR10"] = f"{path}TAIR10_annotation/Araport11_annotation_tidy.gff3"
annotation_files["TAIR10_on_TAIR10"] = f"{path}TAIR10_annotation/TAIR10_annotation_tidy.gff3"

# notice _from_ and _on_ nomenclature
# leave empty if no liftoff or equivalent files
annotation_transfer_files = {}

functional_annotation = {}
functional_annotation["Araport11_on_TAIR10_mapman"] = f"{path}/TAIR10_annotation/functional_annotation/MapMan.Mercator.Araport.gmt"
functional_annotation["Araport11_on_TAIR10_interpro"] = f"{path}/TAIR10_annotation/functional_annotation/Araport11_aegis_on_TAIR10_all_proteins.interpro"
functional_annotation["Araport11_on_TAIR10_eggnog"] = f"{path}/TAIR10_annotation/functional_annotation/out.emapper.annotations"

symbols = {}
symbols["Araport11_on_TAIR10_symbols"] = f"{path}/TAIR10_annotation/Araport11_gene_symbols_def.csv"

genome_dictionary = {}

for tag in genome_files:
    genome_dictionary[tag] = {}