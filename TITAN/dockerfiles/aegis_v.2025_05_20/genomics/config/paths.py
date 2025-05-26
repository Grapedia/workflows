#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os

root = "../../genomes_and_annotation/"


GO_file = f"{root}GO/go.obo"

temp = os.listdir("config/")
all_species = []
for t in temp:
    if t != "paths.py" and ".py" in t:
        all_species.append(t[:-3])