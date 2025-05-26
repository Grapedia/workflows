#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Creating all annotation objects for a single species.

@author: David Navarro
"""

# INPUT: change for other species
from config.PN40024_T2T import genome_files, path

import time
import os
from modules.tools import pickle_load, pickle_save
from modules.genome import Genome
from modules.annotation import Annotation

start = time.time()

aegis_output = f"{path}aegis_output/"
pickle_path = f"{aegis_output}pickles/"
features_path = f"{aegis_output}features/"
gff_path = f"{aegis_output}gffs/"
genome_path = f"{aegis_output}genomes/"
coordinates_path = f"{aegis_output}coordinates/"

os.system(f"mkdir -p {features_path}")
os.system(f"mkdir -p {pickle_path}")
os.system(f"mkdir -p {gff_path}")
os.system(f"mkdir -p {genome_path}")
os.system(f"mkdir -p {coordinates_path}")

genomes = {}

for k, v in genome_files.items():
    genomes[k] = Genome(k, v)

file = f"5.1_on_T2T_ref_annotation.pkl"
a = pickle_load(f"{pickle_path}{file}")

file = f"v3_from_12Xv2_on_T2T_ref_annotation.pkl"
v3 = pickle_load(f"{pickle_path}{file}")

v3.add_gene_names_and_descriptors(file_path=f"../../genomes_and_annotation/grapevine/PN40024/12Xv2_annotation/catalogue/vvi_lab_catalogue_V4_TM39_unofficial_alternative_stats.csv")
v3.remove_noncoding_transcripts_from_coding_genes()
v3.remove_genes_without_names()
v3.id = "catalogue_v3_on_T2T"
v3.name = "catalogue_v3"

a.id = "5.1.1_on_T2T_ref"
a.name = "5.1.1"
a.merge(v3, ignore_overlaps=False, exon_overlap_threshold=0)
a.id = "5.1.1_on_T2T_ref"
a.name = "5.1.1"
a.update_attributes(clean=clean, featurecountsID=True)
a.export_gff(gff_path, UTRs=False)

now = time.time()
lapse = now - start
print(f"\nUpdating to {a.id} annotation took {round((lapse/60)/60, 1)} hours\n")