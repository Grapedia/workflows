#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
2023/06/13

Clean unmet parent references from gff (v3 in particular) and
add gene names from Catalogue

@author: David Navarro
"""

import os

catalogue_path = ("../../listas_de_genes/"
                  "Vv/vvi_lab_catalogue_V3_TM37.csv")

# create folder and add v3 files in there
v3_path = ("../../genomes_and_annotation/grapevine/PN40024/"
           "/T2T_annotation_transfer/"
           "v3_liftoffs_unmapped_parents/")

out_path = ("../../genomes_and_annotation/grapevine/PN40024/"
           "PN40024_helfensteiner/T2T_annotation_transfer/")

v3_files = os.listdir(v3_path)

v3_to_symbol = {}

f_in = open(f"{catalogue_path}", "r", encoding="utf-8")
x = -1
for line in f_in:
    x += 1
    if x > 0:
        temp = line.strip().split("\t")
        if temp[1] not in v3_to_symbol:
            v3_to_symbol[temp[1]] = temp[0]
f_in.close()

for file in v3_files:
    print(f"Reading IDs from {file}")
    f_in = open(f"{v3_path}{file}", "r", encoding="utf-8")
    IDs = []
    for line in f_in:
        if not line.startswith("#"):
            line = line.split("\t")
            attributes = line[-1].split(";")
            for a in attributes:
                if "ID=" in a:
                    ID = a.split("=")[1]
                    if ID not in IDs:
                        IDs.append(ID)
    f_in.close()

    print("Removing unmapped parents and adding gene names to 'mRNA/gene' "
          f"features for {file}")
    f_in = open(f"{v3_path}{file}", "r", encoding="utf-8")
    out = ""
    for line in f_in:
        if line.startswith("#"):
            out += line
            continue
        temp = line.strip().split("\t")
        attributes = temp[-1].split(";")
        parent_found = False
        name_found = False
        ID = ""
        for n, a in enumerate(attributes):
            if a.startswith("ID="):
                ID = a.split("=")[1]
                if "_" in ID:
                    extra_copies = True
                    copy_num = ID.split("_")[1]
                    ID = ID.split("_")[0]
                else:
                    extra_copies = False
                ID = ID.split(".")[0]
            if a.startswith("Parent="):
                temp_a = a.split("=")[1]
                temp_a = temp_a.split(",")
                new_parents = []
                for parent in temp_a:
                    if parent in IDs and parent not in new_parents:
                        new_parents.append(parent)
                parent_index = n
                parents_term = "Parent=" + ",".join(new_parents)
                parent_found = True
            elif "Name=" in a and (temp[2] == "gene" or temp[2] == "mRNA"):
                if ID in v3_to_symbol:
                    if extra_copies:
                        name = f"{v3_to_symbol[ID]}_catalogue_copy{copy_num}"
                    else:
                        name = v3_to_symbol[ID]
                    name_index = n
                    name_found = True
                    names_term = f"Name={name}"   
        if parent_found:
            attributes[parent_index] = parents_term
        if name_found:
            attributes[name_index] = names_term
        elif ID in v3_to_symbol:
            if extra_copies:
                name = f"{v3_to_symbol[ID]}_catalogue_copy{copy_num}"
            else:
                name = v3_to_symbol[ID]
            if temp[2] == "gene":
                attributes.insert(1, f"Name={name}")
            elif temp[2] == "mRNA":
                attributes.insert(2, f"Name={name}")
        temp[-1] = ";".join(attributes)
        out += "\t".join(temp) + "\n"
    f_in.close()

    f_out = open(f"{out_path}{file[:-5]}_clean.gff3", "w", encoding="utf-8")
    f_out.write(out)
    f_out.close()
