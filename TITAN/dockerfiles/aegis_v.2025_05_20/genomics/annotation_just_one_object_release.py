#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Creating all annotation objects for a single species.

@author: David Navarro
"""

# INPUT: change for other species
from config.PN40024 import genome_files, annotation_files, annotation_transfer_files, features_path, pickle_path, species

import time
import os
from modules.tools import pickle_load, pickle_save
from modules.genome import Genome
from modules.annotation import Annotation

start = time.time()
os.system(f"mkdir -p {features_path}")
os.system(f"mkdir -p {pickle_path}")

genomes = {}

for k, v in genome_files.items():
    genomes[k] = Genome(k, v)
    genomes[k].export_feature_sizes()

k = "final_reduced_on_40X_ref"

location = f"{pickle_path}{k}_annotation.pkl"
if os.path.isfile(location):
    print(f"Loading {k} annotation object\n")
    a = pickle_load(location)
else:
    print(f"{k} annotation does not exist")
    annot_name, genome = k.split("_on_")
    a = Annotation(annot_name, annotation_files[k], genomes[genome])
    pickle_save(location, a)

a.generate_sequences(just_CDSs=True, genome=genomes[a.genome])
a.clear_overlaps()
a.detect_gene_overlaps()
a.update_stats(export=True, genome=genomes[a.genome])
a.export_equivalences(stringent=False, return_df=False, export_csv=True, export_self=True, NAs=False)
a.export_gff()

a.combine_transcripts(genomes[a.genome])
a.name = "final_reduced_combined"
a.id = "final_reduced_combined_on_40X_ref"
a.clear_overlaps()
a.detect_gene_overlaps()
a.update_stats(export=True, genome=genomes[a.genome])
a.export_equivalences(stringent=False, return_df=False, export_csv=True, export_self=True, NAs=False)
a.export_gff()

k = "v3_on_12Xv2"
location = f"{pickle_path}{k}_annotation.pkl"
if os.path.isfile(location):
    print(f"Loading {k} annotation object\n")
    a1 = pickle_load(location)
else:
    print(f"{k} annotation does not exist")
    annot_name, genome = k.split("_on_")
    a1 = Annotation(annot_name, annotation_files[k], genomes[genome])
    a1.update_stats(export=True, genome=genomes[a.genome])
    pickle_save(location, a1)

location = f"{pickle_path}catalogue_v3_on_12Xv2_annotation.pkl"
if os.path.isfile(location):
    a2 = pickle_load(location)
else:
    a2 = a1.copy()
    a2.add_gene_names_and_descriptors(file_path=f"../../genomes_and_annotation/grapevine/PN40024/12Xv2_annotation/catalogue/vvi_lab_catalogue_V4_TM39_unofficial_alternative_stats.csv")
    a2.remove_noncoding_transcripts_from_coding_genes()
    a2.remove_genes_without_names()
    a2.id = "catalogue_v3_on_12Xv2"
    a2.name = "catalogue_v3"
    a2.export_gff()
    a2.update_stats(export=True, genome=genomes[a2.genome])
    pickle_save(location, a2)

id_prefix = "Vitvi04_01"
source = "Nemesis"
name = "40X.v4.4"
id = "40X.v4.4_on_40X"

print(a)

a.merge(a2, ignore_overlaps=False, exon_overlap_threshold=0)
a.release(name=name, id=id, source_name=source, id_prefix=id_prefix, spacer=100)

now = time.time()
lapse = now - start
print(f"\nCreating and updating {species} annotation objects as well as exporting all features and gene lengths took {round((lapse/60)/60, 1)} hours\n")