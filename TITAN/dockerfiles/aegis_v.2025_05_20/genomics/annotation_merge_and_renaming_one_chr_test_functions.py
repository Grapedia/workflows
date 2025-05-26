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

# import _pydevd_bundle.pydevd_constants
# _pydevd_bundle.pydevd_constants.PYDEVD_WARN_EVALUATION_TIMEOUT = 30.

start = time.time()
os.system(f"mkdir -p {features_path}")
os.system(f"mkdir -p {pickle_path}")

genomes = {}

for k, v in genome_files.items():
    genomes[k] = Genome(k, v)
    genomes[k].export_feature_sizes()

#"chr01_geneid_on_40X_ref","chr01_glimmer_on_40X_ref", "chr01_transcriptomes_on_40X_ref","chr01_arabidopsis_on_40X_ref", 

annot_tags = ["genemark_on_40X_ref", "augustus_on_40X_ref", "v3_from_12Xv2_raw_on_40X_ref", "psiclass_stranded_on_40X_ref", "psiclass_unstranded_on_40X_ref", "stringtie_default_on_40X_ref", "stringtie_morus_on_40X_ref", "stringtie_stranded_default_on_40X_ref", "stringtie_stranded_morus_on_40X_ref", "stringtie_unstranded_on_40X_ref", "stringtie_unstranded_morus_on_40X_ref"]

annotations = {}

chosen_chromosomes = ["chr01", "chr02", "chr03", "chr04", "chr05", "chr06", "chr07", "chr08", "chr09", "chr10", "chr11", "chr12", "chr13"]

#chosen_chromosomes = ["chr01", "chr02", "chr03", "chr04", "chr05", "chr06", "chr07", "chr08", "chr09", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chrUn"]

for chosen_chromosome in chosen_chromosomes:

    for k, v in annotation_files.items():
        if k in annot_tags:
            location = f"{pickle_path}{chosen_chromosome}_{k}_annotation.pkl"
            if os.path.isfile(location):
                print(f"Loading {k} annotation object\n")
                annotations[k] = pickle_load(location)
            else:
                annot_name, genome = k.split("_on_")
                annotations[k] = Annotation(annot_name, v, genomes[genome], chosen_chromosome=chosen_chromosome)
                annotations[k].export_gff()
                pickle_save(location, annotations[k])

    location = f"{pickle_path}{chosen_chromosome}_abinitio_transcriptome_merge_on_40X_ref_annotation.pkl"
    if os.path.isfile(location):
        b = pickle_load(location)
    else:
        b = annotations[annot_tags[0]].copy()
        del annotations[annot_tags[0]]

        for a in annotations.values():
            b.merge(a, ignore_overlaps=False, exon_overlap_threshold=0)

        b.name = f"{chosen_chromosome}_abinitio_transcriptome_merge"
        b.id = f"{chosen_chromosome}_abinitio_transcriptome_merge_on_40X_ref"

        b.clear_overlaps()
        b.detect_gene_overlaps()
        b.export_equivalences(stringent=False, return_df=False, export_csv=True, export_self=True, NAs=False)

        b.clear_overlaps()
        b.detect_gene_overlaps_old()
        b.name += "_old_overlaps"
        b.id = f"{b.name}_on_40X_ref"
        b.export_equivalences(stringent=False, return_df=False, export_csv=True, export_self=True, NAs=False)
        pickle_save(location, b)


print('End')

now = time.time()
lapse = now - start
print(f"\nCreating and updating {species} annotation objects as well as other procedures took {round((lapse/60)/60, 1)} hours\n")