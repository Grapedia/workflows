#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 26 May 2023

File for Laurent Deluc

@author: David Navarro
"""

path = "../../listas_de_genes/Vv/equivalences/"

v4_file = "all_annotations_compared_to_v4_3_ref.csv"
mapman = "v3_to_v1_description_and_mapman.tsv"

header = "v4_3\tv3\tnon_curated_functional_annotation_transferred_from_v1\tMapMan_bin_from_v3\n"

v4_to_v3 = {}
v3_to_NA = []
v3_to_desc_mapman = {}

f_in = open(path + v4_file, "r", encoding="utf-8")
x = -1
for line in f_in:
    x += 1
    if x > 0:
        line = line.strip().split("\t")
        origin = line[3]
        v4 = line[0]
        query = line[1]
        if v4 != "NA":
            if v4 not in v4_to_v3:
                if query != "NA" and origin == "v3_liftoff":                
                    v4_to_v3[v4] = [query]
                else:
                    v4_to_v3[v4] = []
            elif origin == "v3_liftoff" and query != "NA":
                v4_to_v3[v4].append(query)
        elif origin == "v3_liftoff" and query != "NA":
            v3_to_NA.append(query)
f_in.close()

f_in = open(path + mapman, "r", encoding="utf-8")
x = -1
for line in f_in:
    x += 1
    if x > 0:
        line = line.strip().split("\t")
        v3 = line[0]
        desc = line[1]
        mpm = line[2]
        if v3 not in v3_to_desc_mapman:
            v3_to_desc_mapman[v3] = f"{desc}\t{mpm}"
        else:
            print(f"{v3} not unique")
f_in.close()

text_out = header
included_v3s = []

for v4 in v4_to_v3:
    if v4_to_v3[v4] != []:
        for v3 in v4_to_v3[v4]:
            if v3 in v3_to_desc_mapman:
                desc_mapman = v3_to_desc_mapman[v3]
            else:
                desc_mapman = "NA\tNA"
            text_out += f"{v4}\t{v3}\t{desc_mapman}\n"
            included_v3s.append(v3)
    else:
        text_out += f"{v4}\tNA\tNA\tNA\n"

for v3 in v3_to_NA:
    if v3 in v3_to_desc_mapman:
        desc_mapman = v3_to_desc_mapman[v3]
    else:
        desc_mapman = "NA\tNA"
    text_out += f"NA\t{v3}\t{desc_mapman}\n"
    included_v3s.append(v3)


for v3 in v3_to_desc_mapman:
    if v3 not in included_v3s:
        desc_mapman = v3_to_desc_mapman[v3]
        text_out += f"NA\t{v3}\t{desc_mapman}\n"

f_out = open(path + "v4_3_to_v3_to_description_and_mapman.tsv", "w", 
             encoding="utf-8")
f_out.write(text_out)
f_out.close()

