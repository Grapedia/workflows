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

k = "5.1_marked_on_T2T_ref_minus_non_coding_aegis_dapfit"

annotation_gff = f"{gff_path}{k}.gff3"

location = f"{pickle_path}{k}_annotation.pkl"
if os.path.isfile(location):
    a = pickle_load(location)
else:
    annot_name = k.split("_on_")[0]
    genome = "T2T_ref"
    a = Annotation(annot_name, annotation_gff, genomes[genome])
    pickle_save(location, a)

TE_file = "../../genomes_and_annotation/grapevine/PN40024/T2T_annotation/transposons/interpro_TEs.csv"

a.export_gff(gff_path2, UTRs=False)
a.mark_transposable_element_genes(TE_file)
a.remove_genes_with_small_CDSs()
a.remove_TE_genes()
a.update_attributes(clean=True)
a.update_stats(custom_path=stats_path, export=True, genome=genomes[a.genome], max_x=1000)
a.export_gff(gff_path2, UTRs=False)

now = time.time()
lapse = now - start
print(f"\nUpdating to {a.id} annotation took {round((lapse/60)/60, 1)} hours\n")