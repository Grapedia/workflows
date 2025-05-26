#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Creating all annotation objects for a single species.

@author: David Navarro
"""

# INPUT: change for other species
from config.PN40024 import genome_files, annotation_files, path, species

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
os.system(f"mkdir -p {pickle_path}")
os.system(f"mkdir -p {gff_path}")
os.system(f"mkdir -p {features_path}")

genomes = {}

for k, v in genome_files.items():
    genomes[k] = Genome(k, v)
    genomes[k].export_feature_sizes()

k = "5.1_on_T2T_ref"

location = f"{pickle_path}{k}_annotation.pkl"
if os.path.isfile(location):
    print(f"Loading {k} annotation object\n")
    a = pickle_load(location)
else:
    print(f"{k} annotation does not exist")
    annot_name, genome = k.split("_on_")
    a = Annotation(annot_name, annotation_files[k], genomes[genome])
    pickle_save(location, a)

a.add_gene_symbols("../../genomes_and_annotation/grapevine/PN40024/T2T_annotation/functional_annotation/grapevine_catalogue_v3.xlsx")
a.remove_genes_without_symbols()
a.update_attributes(symbols=False, symbols_as_descriptors=True)
a.export_gff(custom_path=gff_path, UTRs=True)

now = time.time()
lapse = now - start
print(f"\nCreating and updating {species} annotation objects as well as exporting all features and gene lengths took {round((lapse/60)/60, 1)} hours\n")