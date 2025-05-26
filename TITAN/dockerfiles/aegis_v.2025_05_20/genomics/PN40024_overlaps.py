#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 28 2023

Reading annotation objects, calculating and exporting overlaps.

@author: David Navarro
"""

import time
import os
from modules.tools import pickle_load
from modules.genefunctions import export_group_equivalences
from modules.geneclasses import Genome, Annotation

# INPUT
start = time.time()
root = "../../genomes_and_annotation/grapevine/PN40024/"
transfer_path = f"{root}40X_annotation_transfer/"
pickle_path = f"{root}pickle_objects/"
overlaps_path = f"{root}overlaps/"

# PROGRAM
ref_annotation_tags = ["v4.3_on_40X_ref","8x_on_40X_ref", "NCBI_on_40X_ref",
                       "v0_on_40X_ref", "v1_on_40X_ref", "v2_on_40X_ref",
                       "v3_on_40X_ref"]

alt_annotation_tags = ["v4.1_on_40X_alt", "8x_on_40X_alt", 
                       "NCBI_on_40X_alt", "v0_on_40X_alt", "v1_on_40X_alt",
                       "v2_on_40X_alt", "v3_on_40X_alt"]

annotations = {}

for n, tag in enumerate(ref_annotation_tags):
    location = f"{pickle_path}{tag}_annotation.pkl"
    if os.path.isfile(location):
        print(f"Loading {ref_annotation_tags[n]} annotation object\n")
        annotations[tag] = pickle_load(location)

# multidirectional attempts between all annotations on REF (also does single)
# export_group_equivalences(list(annotations.values()), overlaps_path, True, 
#                           clear_overlaps=False)
export_group_equivalences(list(annotations.values()), overlaps_path, True,
                          stringent=False, clear_overlaps=False)

annotations.clear()

for n, tag in enumerate(alt_annotation_tags):
    location = f"{pickle_path}{tag}_annotation.pkl"
    if os.path.isfile(location):
        print(f"Loading {alt_annotation_tags[n]} annotation object\n")
        annotations[tag] = pickle_load(location)

# non multidirectional (only to target) using 4.1 ALT
# annotations["v4.1_on_40X_alt"].target = True
# export_group_equivalences(list(annotations.values()), clear_overlaps=False)

# # multidirectional attempts between all annotaions on ALT
# export_group_equivalences(list(annotations.values()), overlaps_path, True,
#                           clear_overlaps=False)
export_group_equivalences(list(annotations.values()), overlaps_path, True,
                          stringent=False, clear_overlaps=False)

now = time.time()
lapse = now - start
print("\nGenerating and exporting all individual and group PN40024 annotation "
      f"overlaps took {round((lapse/60)/60, 1)} hours\n")