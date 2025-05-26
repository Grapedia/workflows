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
from modules.annotation import Annotation

from config.PN40024 import *

start = time.time()

aegis_output = f"{path}aegis_output/"
pickle_path = f"{aegis_output}pickles/"
overlaps_path = f"{aegis_output}overlaps/"

os.system(f"mkdir -p {overlaps_path}")

annotation_files.update(annotation_transfer_files)
annotations = {}

annotations_to_process = ['5.0_on_T2T_ref', '5.1_on_T2T_ref', 'v4.3_from_40X_ref_on_T2T_ref', 'v3_from_12Xv2_on_T2T_ref', '8x_from_8x_on_T2T_ref', '12XNCBI_from_12XNCBI_on_T2T_ref', 'v0_from_12Xv2_on_T2T_ref', 'v1_from_12Xv2_on_T2T_ref', 'v2_from_12Xv2_on_T2T_ref']

for n, tag in enumerate(annotation_files):
    location = f"{pickle_path}{tag}_annotation.pkl"
    if tag in annotations_to_process:
        if os.path.isfile(location):
            print(f"Loading {tag} annotation object\n")
            annotations[tag] = pickle_load(location)

#multidirectional attempts between all annotations on REF (also does single)
export_group_equivalences(list(annotations.values()), overlaps_path, multidirectional=True, stringent=True, clear_overlaps=False)

now = time.time()
lapse = now - start
print(f"\nGenerating and exporting all individual and group PN40024 T2T annotation overlaps took {round((lapse/60)/60, 1)} hours\n")