#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: David Navarro
"""

root = f"../../genomes_and_annotation/grapevine/PN40024/"

functional_annotation_f = f"{root}functional_annotation_output/5.1_on_T2T_ref_functional_annotation_summary.csv"

gene_coordinates_f = f"{root}T2T_annotation/gene_coordinates/5.1_on_T2T_ref_gene_coordinates.tsv"

out_file = f"{root}functional_annotation_output/5.1_on_T2T_ref_functional_annotation_summary_with_coordinates.tsv"


gene_coordinates = {}

f_in = open(gene_coordinates_f)

for n, line in enumerate(f_in):
    if n > 0:
        line = line.strip().split("\t")
        if line[0] not in gene_coordinates:
            gene_coordinates[line[0]] = (line[1], line[2], line[3])
        else:
            print("Error")

f_in.close()

out = []

f_in = open(functional_annotation_f)

for n, line in enumerate(f_in):
    line = line.strip().split("\t")
    if n == 0:
        line.insert(1, "gene_end")
        line.insert(1, "gene_start")
        line.insert(1, "chromosome")
    else:
        gene = line[0]
        if gene in gene_coordinates:
            line.insert(1, gene_coordinates[gene][2])
            line.insert(1, gene_coordinates[gene][1])
            line.insert(1, gene_coordinates[gene][0])
        else:
            print("Error 2")
    out.append(line)

f_in.close()


out_string = ""

for line in out:
    out_string += "\t".join(line)
    out_string += "\n"

f_out = open(out_file, "w", encoding="utf-8")
f_out.write(out_string)

f_out.close()