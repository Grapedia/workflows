#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
2023/05/22

Changing the incorrect gtf/gff3 format of 8x grapevine annotation into gff3

@author: David Navarro
"""

path_in = ("../../genomes_and_annotation/grapevine/PN40024/"
           "8X_annotation_genoscope/original/")

f_in = open(path_in + "8x_genoscope.gtf", "r", encoding="utf-8")

# there is only 1 mRNA per gene

out = ""
count = 0
for line in f_in:
    temp = line.strip().split("\t")
    ft = temp[2]
    attributes = temp[-1].split(" ; ")
    ID = attributes[0].split(" ")[1]

    if ft == "gene":
        gene_id = ID
        count = 0
        if "Complete 0" in line:
            complete = 0
        elif "Complete 1" in line:
            complete = 1
        else:
            print("code error")
        temp[-1] = f"ID={ID};Complete={complete}"


    if ft == "mRNA":
        if "Complete 0" in line:
            complete = 0
        elif "Complete 1" in line:
            complete = 1
        else:
            print("code error")
        temp[-1] = f"ID={ID};Parent={gene_id};Complete={complete}"

    if ft == "CDS":
        temp[-1] = f"ID={ID}.CDS1;Parent={ID}"

    if ft == "UTR":
        count += 1
        temp[-1] = f"ID={ID}.UTR{count};Parent={ID}"


    out += "\t".join(temp) + "\n"
f_in.close()


f_out = open(("/home/tomslab3/Escritorio/postdoc_i2sysbio/"
              "genomes_and_annotation/PN40024_helfensteiner/"
              "8X_annotation_genoscope/original/8x_corrected.gff3"),
              "w", encoding="utf-8")
f_out.write(out)
f_out.close()

