#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Creating all annotation objects for a single species.

@author: David Navarro
"""

# INPUT: change for other species
from config.PN40024 import genome_files, path, annotation_files

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
gff_path2 = f"{aegis_output}gffs_test/"
genome_path = f"{aegis_output}genomes/"
coordinates_path = f"{aegis_output}coordinates/"
stats_path = f"{aegis_output}stats/"

os.system(f"mkdir -p {features_path}")
os.system(f"mkdir -p {pickle_path}")
os.system(f"mkdir -p {gff_path}")
os.system(f"mkdir -p {genome_path}")
os.system(f"mkdir -p {coordinates_path}")
os.system(f"mkdir -p {stats_path}")
os.system(f"mkdir -p {gff_path2}")


genomes = {}

for k, v in genome_files.items():
    genomes[k] = Genome(k, v)

k = "5.1_on_T2T_ref"

annotation_gff = f"{gff_path}{k}.gff3"

location = f"{pickle_path}{k}_annotation.pkl"
if os.path.isfile(location):
    a = pickle_load(location)
else:
    annot_name = k.split("_on_")[0]
    genome = "T2T_ref"
    a = Annotation(annot_name, annotation_gff, genomes[genome])
    pickle_save(location, a)


print(a.id, a.suffix, a.feature_suffix)
a.export_gff(custom_path=gff_path, UTRs=False)

nc = a.copy()
nc.remove_coding_genes_and_transcripts()
print(nc.id, nc.suffix, nc.feature_suffix)
nc.export_gff(custom_path=gff_path, UTRs=False)

a.remove_non_coding_genes_and_transcripts()
print(a.id, a.suffix, a.feature_suffix)
a.export_gff(custom_path=gff_path, UTRs=False)
a.export_all_features(feature_output="all", promoters=False, verbose=False, path=features_path)

now = time.time()
lapse = now - start
print(f"\nUpdating to {a.id} annotation took {round((lapse/60)/60, 1)} hours\n")