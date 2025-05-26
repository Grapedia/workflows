#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
2023/05/22

Changing the chromosome numbers of a fasta file

@author: David Navarro
"""

path_in = ("../../genomes_and_annotation/grapevine/PN40024/"
           "8X_genome_genoscope/original/")

f_in = open(path_in + "VV_chr8x.fsa", "r", encoding="utf-8")

out = ""
for line in f_in:
    if line[0] == ">":
        temp = line.split(" ")
        if len(temp) > 1:
            if "chr" in temp[-1]:
                line = f">{temp[-1]}"
            else:
                print("code error")
    out += line
f_in.close()


f_out = open(("/home/tomslab3/Escritorio/postdoc_i2sysbio/"
              "genomes_and_annotation/PN40024_helfensteiner/"
              "8X_genome_genoscope/8x.fasta"), "w", encoding="utf-8")
f_out.write(out)
f_out.close()

