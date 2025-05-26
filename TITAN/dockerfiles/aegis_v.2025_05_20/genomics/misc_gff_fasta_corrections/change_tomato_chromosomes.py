#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
2023/05/22

Changing the chromosome numbers of a gff.

@author: David Navarro
"""

path= "../../genomes_and_annotation/tomato/Heinz_ITAG4_annotation/"

f_in = open(f"{path}ITAG4.0_priority_4.2_merge_tidy_fc_myb24_mod.gff3", "r",
            encoding="utf-8")
out = ""
for line in f_in:
    if line[0] != "#":
        line = line.split("\t")
        line[0] = "chr" + line[0].split("h")[1]
        line = "\t".join(line)
    out += line
f_in.close()

f_out = open(f"{path}ITAG4.0_priority_4.2_merge_tidy_fc_myb24_mod_DAP.gff3", 
             "w", encoding="utf-8")
f_out.write(out)
f_out.close()
