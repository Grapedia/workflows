#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
2023/06/04

v4_1 IDs

@author: David Navarro
"""

v4_1_alt = ("../../genomes_and_annotation/grapevine/PN40024/40X_annotation/"
            "original/v4_1_alt.gff3")

non_alt_IDs_out = ("../../genomes_and_annotation/grapevine/PN40024/"
                   "40X_annotation/original/v4_1_alt_non_alt_IDs.txt")

v4_3_1_merge = ("../../genomes_and_annotation/grapevine/PN40024/"
                "40X_annotation/incorrect_merge/v4_3_ref_1_alt_merge.gff3")

duplicated_merge_IDs_out = ("../../genomes_and_annotation/grapevine/PN40024"
                            "/40X_annotation/incorrect_merge/"
                            "v4_3_ref_1_alt_merge_duplicated_IDs.txt")

v4_1_merge = ("../../genomes_and_annotation/grapevine/PN40024/"
              "40X_annotation/original/v4_1_merged.gff3")

duplicated_4_1_merge_IDs_out = ("../../genomes_and_annotation/grapevine/"
                                "PN40024/40X_annotation/original/"
                                "v4_1_merged_duplicated_IDs.txt")


non_alt_IDs = []

f_in = open(v4_1_alt, "r", encoding="utf-8")
for line in f_in:
    if line[0] != "#":
        line = line.strip().split("\t")
        ft = line[2]
        if ft == "gene" or ft == "pseudogene":
            attributes = line[-1].split(";")
            for a in attributes:
                if "ID=" in a:
                    if "_alt" not in a:
                        ID = a.split("=")[1]
                        if ID not in non_alt_IDs:
                            non_alt_IDs.append(ID)

f_in.close()

print(f"There are {len(non_alt_IDs)} non-alt IDs in v4_1_alt")

non_alt_IDs.sort()
f_out = open(non_alt_IDs_out, "w", encoding="utf-8")
f_out.write("\n".join(non_alt_IDs))
f_out.close()


IDs = []
duplicated_merge_IDs = []

f_in = open(v4_3_1_merge, "r", encoding="utf-8")
for line in f_in:
    if line[0] != "#":
        line = line.strip().split("\t")
        ft = line[2]
        if ft == "gene" or ft == "pseudogene":
            attributes = line[-1].split(";")
            for a in attributes:
                if "ID=" in a:
                    ID = a.split("=")[1]
                    if ID in IDs:
                        duplicated_merge_IDs.append(ID)
                    IDs.append(ID)
f_in.close()

print(f"There are {len(duplicated_merge_IDs)} duplicated IDs in "
      "v4_3_ref_1_alt_merge")
duplicated_merge_IDs.sort()
f_out = open(duplicated_merge_IDs_out, "w", encoding="utf-8")
f_out.write("\n".join(duplicated_merge_IDs))
f_out.close()

IDs = []
duplicated_merge_IDs = []

f_in = open(v4_1_merge, "r", encoding="utf-8")
for line in f_in:
    if line[0] != "#":
        line = line.strip().split("\t")
        ft = line[2]
        if ft == "gene" or ft == "pseudogene":
            attributes = line[-1].split(";")
            for a in attributes:
                if "ID=" in a:
                    ID = a.split("=")[1]
                    if ID in IDs:
                        duplicated_merge_IDs.append(ID)
                    IDs.append(ID)
f_in.close()

print(f"There are {len(duplicated_merge_IDs)} duplicated IDs in "
      "v4_1_merged")
duplicated_merge_IDs.sort()

if len(duplicated_4_1_merge_IDs_out) > 0:
    f_out = open(duplicated_4_1_merge_IDs_out, "w", encoding="utf-8")
    f_out.write("\n".join(duplicated_merge_IDs))
    f_out.close()