#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Renaming CYP ids to check the issue.

@author: David Navarro
"""

# INPUT: change for other species
from config.PN40024 import *

import time
import os
from modules.tools import pickle_load, pickle_save
from modules.annotation import Annotation


start = time.time()
aegis_output = f"{path}aegis_output/"
pickle_path = f"{aegis_output}pickles/"
features_path = f"{aegis_output}features/"
gff_path = f"{aegis_output}gffs/"
genome_path = f"{aegis_output}genomes/"
coordinates_path = f"{aegis_output}coordinates/"



k = "cyps_on_12Xv2"

location = f"{pickle_path}{k}_annotation.pkl"
if os.path.isfile(location):
    print(f"Loading {k} annotation object\n")
    a = pickle_load(location)


a.rename_ids(prefix="cyp", spacer=1, correspondences=True, custom_path=gff_path)
a.id = "cyps_renamed_on_12Xv2"
a.name = "cyps_renamed"
a.export_gff(skip_atypical_fts=True, UTRs=False, custom_path=gff_path)

now = time.time()
lapse = now - start
print(f"\nTook {round((lapse/60)/60, 1)} hours\n")