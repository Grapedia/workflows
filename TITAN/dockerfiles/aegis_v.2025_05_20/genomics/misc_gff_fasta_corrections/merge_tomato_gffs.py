#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
2023/05/19

Merging 4.0 and 4.1 annotations from Tomato in two ways

@author: David Navarro
"""

from re import search

file_4_0 = "../../genomes_and_annotation/tomato/Heinz_ITAG4_annotation/original/ITAG4.0_tidy.gff3"
file_4_1 = "../../genomes_and_annotation/tomato/Heinz_ITAG4_annotation/original/ITAG4.1_tidy.gff3"
file_4_2 = "../../genomes_and_annotation/tomato/Heinz_ITAG4_annotation/original/ITAG4.2_beta_tidy.gff3"
file_4_0_out = "../../genomes_and_annotation/tomato/Heinz_ITAG4_annotation/original/ITAG4.0_priority_4.1_merge_v2.gff3"
file_4_1_out = "../../genomes_and_annotation/tomato/Heinz_ITAG4_annotation/original/ITAG4.1_priority_4.0_merge_v2.gff3"
# SINCE V4.1 added no genes, we renamed the file
file_4_2_out = "../../genomes_and_annotation/tomato/Heinz_ITAG4_annotation/original/ITAG4.0_priority_4.2_merge.gff3"

d_4_0 = {}
d_4_1 = {}
d_4_2 = {}

f_in = open(file_4_0, "r", encoding="utf-8")
for line in f_in:
    if line.startswith("#"):
        continue
    m = search(r"Solyc[0-9]{2}g[0-9]{6}", line)
    ID = m.group(0)
    terms = line.strip().split("\t")
    if ID not in d_4_0:
        d_4_0[ID] = {}
        d_4_0[ID]["lines"] = line
    else:
        d_4_0[ID]["lines"] += line
    if terms[2] == "gene":
        d_4_0[ID]["ch"] = terms[0]     
        d_4_0[ID]["coord"] = [terms[4], terms[5]]
f_in.close()


f_in = open(file_4_1, "r", encoding="utf-8")
for line in f_in:
    if line.startswith("#"):
        continue
    m = search(r"Solyc[0-9]{2}g[0-9]{6}", line)
    ID = m.group(0)
    terms = line.strip().split("\t")
    if ID not in d_4_1:
        d_4_1[ID] = {}
        d_4_1[ID]["lines"] = line
    else:
        d_4_1[ID]["lines"] += line
    if terms[2] == "gene":
        d_4_1[ID]["ch"] = terms[0]     
        d_4_1[ID]["coord"] = [terms[4], terms[5]]
f_in.close()

f_in = open(file_4_2, "r", encoding="utf-8")
for line in f_in:
    if line.startswith("#"):
        continue
    m = search(r"Solyc[0-9]{2}g[0-9]{6}", line)
    ID = m.group(0)
    terms = line.strip().split("\t")
    if ID not in d_4_2:
        d_4_2[ID] = {}
        d_4_2[ID]["lines"] = line
    else:
        d_4_2[ID]["lines"] += line
    if terms[2] == "gene":
        d_4_2[ID]["ch"] = terms[0]     
        d_4_2[ID]["coord"] = [terms[4], terms[5]]
f_in.close()


f_in = open(file_4_0, "r", encoding="utf-8")
output_4_0_priority = f_in.read()
f_in.close()

f_in = open(file_4_1, "r", encoding="utf-8")
output_4_1_priority = f_in.read()
f_in.close()


for gene, attributes in d_4_0.items():
    if gene not in d_4_1:
        equivalent_gene = False
        for g in d_4_1.values():
            if (g["ch"] == attributes["ch"] 
                and g["coord"] == attributes["coord"]):
                equivalent_gene = True
                break
        if not equivalent_gene:
            output_4_1_priority += attributes["lines"]

for gene, attributes in d_4_1.items():
    if gene not in d_4_0:
        equivalent_gene = False
        for g in d_4_0.values():
            if (g["ch"] == attributes["ch"] 
                and g["coord"] == attributes["coord"]):
                equivalent_gene = True
                break
        if not equivalent_gene:
            output_4_0_priority += attributes["lines"]

f_out = open(file_4_0_out, "w", encoding="utf-8")
f_out.write(output_4_0_priority)
f_out.close()

f_out = open(file_4_1_out, "w", encoding="utf-8")
f_out.write(output_4_1_priority)
f_out.close()


f_in = open(file_4_0, "r", encoding="utf-8")
output_4_0_priority_also_4_2 = f_in.read()
f_in.close()

print("Genes added from 4.2 genes:")
for gene, attributes in d_4_2.items():
    if gene not in d_4_0:
        equivalent_gene = False
        for g in d_4_0.values():
            if (g["ch"] == attributes["ch"] 
                and g["coord"] == attributes["coord"]):
                equivalent_gene = True
                break
        if not equivalent_gene:
            print(gene)
            output_4_0_priority_also_4_2 += attributes["lines"]        

# NO GENES WERE ADDED from 4.1 that's why we renamed the output file
print("Genes added from 4.1 genes:")
for gene, attributes in d_4_1.items():
    if gene not in d_4_2:
        equivalent_gene = False
        for g in d_4_2.values():
            if (g["ch"] == attributes["ch"] 
                and g["coord"] == attributes["coord"]):
                equivalent_gene = True
                break
        if not equivalent_gene:
            if gene not in d_4_0:
                equivalent_gene = False
                for g in d_4_0.values():
                    if (g["ch"] == attributes["ch"] 
                        and g["coord"] == attributes["coord"]):
                        equivalent_gene = True
                        break
                if not equivalent_gene:
                    print(gene)
                    output_4_0_priority_also_4_2 += attributes["lines"]   


f_out = open(file_4_2_out, "w", encoding="utf-8")
f_out.write(output_4_0_priority_also_4_2)
f_out.close()