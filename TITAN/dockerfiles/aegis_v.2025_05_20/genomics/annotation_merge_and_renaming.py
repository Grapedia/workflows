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
from modules.geneclasses import Genome, Annotation

start = time.time()
os.system(f"mkdir -p {features_path}")
os.system(f"mkdir -p {pickle_path}")

genomes = {}

for k, v in genome_files.items():
    genomes[k] = Genome(k, v)
    genomes[k].export_feature_sizes()

annot_tags = ["transcriptomes_on_40X_ref", "stringtie_morus_on_40X_ref", "stringtie_unstranded_on_40X_ref", "geneid_on_40X_ref", "glimmer_on_40X_ref", "augustus_on_40X_ref", "arabidopsis_on_40X_ref", "exonerate_on_40X_ref", "v3_from_12Xv2_raw_on_40X_ref"]

annotations = {}

for k, v in annotation_files.items():
    if k in annot_tags:
        location = f"{pickle_path}{k}_annotation.pkl"
        if os.path.isfile(location):
            print(f"Loading {k} annotation object\n")
            annotations[k] = pickle_load(location)
        else:
            annot_name, genome = k.split("_on_")
            annotations[k] = Annotation(annot_name, v, genomes[genome])
            pickle_save(location, annotations[k])
        annotations[k].export_gff(UTRs=True)

location = f"{pickle_path}abinitio_transcriptome_merge_on_40X_ref_annotation.pkl"
if os.path.isfile(location):
    b = pickle_load(location)
else:
    b = annotations["transcriptomes_on_40X_ref"].copy()
    del annotations["transcriptomes_on_40X_ref"]

    for a in annotations.values():
        b.merge_annotations(a, overlap_threshold=0)

    b.name = "abinitio_transcriptome_merge"
    b.id = "abinitio_transcriptome_merge_on_40X_ref"
    b.export_gff()
    pickle_save(location, b)

location = f"{pickle_path}full_abinitio_transcriptome_merge_on_40X_ref_annotation.pkl"
if os.path.isfile(location):
    b = pickle_load(location)
else:
    b.make_alternative_transcripts_into_genes()
    b.name = "full_abinitio_transcriptome_merge"
    b.id = "full_abinitio_transcriptome_merge_on_40X_ref"
    b.update_stats(export=True, genome=genomes[b.genome])
    b.export_gff()
    b.export_unique_proteins(genome=genomes[b.genome])
    b.generate_sequences(genome=genomes[b.genome])
    b.detect_gene_overlaps()
    # b.add_blast_hits("Araport", f"{b.path}/blast_results/David_test3_Araport11.diamond")
    # b.add_blast_hits("Viridiplantae", f"{b.path}/blast_results/David_test3_Viridiplantae.diamond")
    # b.add_blast_hits("Eudicotyledons", f"{b.path}/blast_results/David_test3_Eudicotyledons.diamond")
    pickle_save(location, b)

# location = f"{pickle_path}abinitio_transcriptome_merge_reduced_on_40X_ref_annotation.pkl"
# if os.path.isfile(location):
#     b = pickle_load(location)
# else:
#     b.remove_redundancy2(source_priority=["Araport", "Viridiplantae", "Eudicotyledons"])
#     b.name = "abinitio_transcriptome_merge_reduced"
#     b.id = "abinitio_transcriptome_merge_reduced_on_40X_ref"
#     b.export_gff()
#     b.update_stats(export=True, genome=genomes[b.genome])
#     b.detect_gene_overlaps()
#     b.export_equivalences(stringent=False, return_df=False, export_csv=True, export_self=True, NAs=False)
#     pickle_save(location, b)

now = time.time()
lapse = now - start
print(f"\nCreating and updating {species} annotation objects as well as "
      f"other procedures took {round((lapse/60)/60, 1)} hours\n")