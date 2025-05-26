#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from config.paths import root, GO_file

species="grapevine"
# assembly/annotation folder if any
assembly="malbec/"
species_path=f"{root}{species}/"
path = f"{species_path}{assembly}"

genome_files = {}
genome_files["malbec_magde"] = f"{path}magdeleine_assembly/Vitis_malbecmagdeleine.Malbec_Magdeleine_v1.0.dna.toplevel.fa"
genome_files["malbec_prune"] = f"{path}prunelard_assembly/Vitis_malbecprunelard.Malbec_Prunelard_v1.0.dna.toplevel.fa"
genome_files["malbec_prune_changed"] = f"{path}prunelard_assembly/Malbec_prunelard_changed.fa"

annotation_files = {}
annotation_files["v1_on_malbec_magde"] = f"{path}magdeleine_assembly/vitis_malbecmagdeleine.gff"
annotation_files["PN40024.v3_on_malbec_magde"] = f"{path}magdeleine_assembly/liftoff/Malbec_magde_v3.gff3"

annotation_files["v1_on_malbec_prune"] = f"{path}prunelard_assembly/vitis_malbecprunelard.gff"
annotation_files["v1_on_malbec_prune_changed"] = f"{path}prunelard_assembly/vitis_malbecprunelard_changed.gff"
#annotation_files["PN40024.v3_on_malbec_prune"] = f"{path}prunelard_assembly/liftoff/Malbec_prune_v3.gff3"
annotation_files["PN40024.v3_on_malbec_prune_changed"] = f"{path}prunelard_assembly/liftoff/Malbec_prune_v3_changed.gff3"
annotation_files["titan_on_malbec_prune_changed"] = f"{path}prunelard_assembly/final_annotation_testMalbec_1.gff3"

genome_dictionary = {}

for tag in genome_files:
    genome_dictionary[tag] = {}