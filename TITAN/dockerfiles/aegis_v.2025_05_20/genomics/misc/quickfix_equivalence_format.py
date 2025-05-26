#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 4 Jul 2023

Fixing issues with equivalences format

@author: David Navarro
"""

import os
import pandas as pd
import io

out_path = "../../genomes_and_annotation/grapevine/PN40024/overlaps/old/format_corrected/"
original_path = "../../genomes_and_annotation/grapevine/PN40024/overlaps/old/format_corrected/original/"
v4_path = "../../genomes_and_annotation/grapevine/PN40024/overlaps/old/format_corrected/v4/"
original_files = os.listdir(original_path)
v4_files = os.listdir(v4_path)
final_path = "../../genomes_and_annotation/grapevine/PN40024/overlaps/old/format_corrected/final/"


for file in original_files:
    if ".csv" in file and not file.startswith("."):
        out = ""
        out_file = file.replace("_queryNAs_verbose_targetNAs", "")
        out_file = out_file.replace("other_", "")
        out_file = out_file.replace("liftoff_alt", "40X_alt")
        out_file = out_file.replace("liftoff_ref", "40X_ref")
        f_in = open(original_path + file, "r", encoding="utf-8")
        x = 0
        for line in f_in:
            x += 1
            line = line.strip().split("\t")
            if x == 1:
                line[0] = line[0].replace("_liftoff_alt", "")
                line[0] = line[0].replace("_liftoff_ref", "")
                target_origin = line[0]
                line[0] = "target_id"
                line[2] = "target_synteny_conserved"
                line[3] = "query_synteny_conserved"
                line[4] = "same_strand"
                line[-2] = "query_origin"
                line[-1] = "overlap_score"
                line.insert(-2, "target_origin")
            else:
                swap = line[3]
                line[3] = line[2]
                line[2] = swap
                line[-2] = line[-2].replace("_liftoff_alt", "")
                line[-2] = line[-2].replace("_liftoff_ref", "")
                if (line[-1] == "NA" and "unmapped" not in line[-2] 
                                     and "unampped" not in line[-2]):
                    line[-1] = "0"
                line[-2] = line[-2].replace("_unmapped", "")
                line[-2] = line[-2].replace("_unampped", "")
                line[-2] = line[-2].replace("v4_1_alt", "v4.1")
                line[-2] = line[-2].replace("v4_3_ref", "v4.3")
                if line[2] == "None":
                    line[2] = "NA"
                if line[3] == "None":
                    line[3] = "NA"
                if line[5] != "NA":
                    line[5] = str(round(float(line[5]), 1))
                if line[7] != "NA":
                    line[7] = str(round(float(line[7]), 1))
                if line[9] != "NA":
                    line[9] = str(round(float(line[9]), 1))
                if line[-1] != "NA":
                    line[-1] = str(int(line[-1]))
                line.insert(-2, target_origin)
                if line[0] == "NA":
                    line[-3] = "NA"
                if line[1] == "NA":
                    line[-2] = "NA"
            del line[8]
            del line[6]
            line = "\t".join(line) + "\n"
            out += line
        f_in.close()
        f_out = open(out_path + out_file, "w", encoding="utf-8")
        f_out.write(out)
        f_out.close()

for file in v4_files:
    if ".csv" in file and not file.startswith("."):
        out = ""
        out_file = file.replace("_queryNAs_verbose_targetNAs", "")
        out_file = out_file.replace("other_", "")
        out_file = out_file.replace("alt", "40X_alt")
        out_file = out_file.replace("ref", "40X_ref")
        f_in = open(v4_path + file, "r", encoding="utf-8")
        x = 0
        for line in f_in:
            x += 1
            line = line.strip().split("\t")
            if x == 1:
                line[0] = line[0].replace("v4_1_alt", "v4.1")
                line[0] = line[0].replace("v4_3_ref", "v4.3")
                target_origin = line[0]
                line[0] = "target_id"
                line.insert(2, "target_synteny_conserved")
                line[3] = "query_synteny_conserved"
                line[4] = "same_strand"
                line[-2] = "query_origin"
                line[-1] = "overlap_score"
                line.insert(-2, "target_origin")
            else:
                line.insert(2, "NA")
                line[-2] = line[-2].replace("_liftoff_alt", "")
                line[-2] = line[-2].replace("_liftoff_ref", "")
                if (line[-1] == "NA" and "unmapped" not in line[-2] 
                                     and "unampped" not in line[-2]):
                    line[-1] = "0"
                line[-2] = line[-2].replace("_unmapped", "")
                line[-2] = line[-2].replace("_unampped", "")
                line[-2] = line[-2].replace("v4_1_alt", "v4.1")
                line[-2] = line[-2].replace("v4_3_ref", "v4.3")
                if line[2] == "None":
                    line[2] = "NA"
                if line[3] == "None":
                    line[3] = "NA"
                if line[5] != "NA":
                    line[5] = str(round(float(line[5]), 1))
                if line[7] != "NA":
                    line[7] = str(round(float(line[7]), 1))
                if line[9] != "NA":
                    line[9] = str(round(float(line[9]), 1))
                if line[-1] != "NA":
                    line[-1] = str(int(line[-1]))                                
                line.insert(-2, target_origin)
                if line[0] == "NA":
                    line[-3] = "NA"
                if line[1] == "NA":
                    line[-2] = "NA"
            del line[8]
            del line[6]
            line = "\t".join(line) + "\n"
            out += line
        f_in.close()
        f_out = open(out_path + out_file, "w", encoding="utf-8")
        f_out.write(out)
        f_out.close()

corrected_files = os.listdir(out_path)

out = ""
x = 0
for file in corrected_files:
    if "ref" in file and "filtered" in file:
        x += 1
        out_file = "all_40X_ref_equivalences_filtered.csv"
        j = 0
        f_in = open(out_path + file, "r", encoding="utf-8")
        for line in f_in:
            append = False
            line = line.strip().split("\t")
            j += 1
            if j == 1:
                if x == 1:
                    append = True
            elif line[0] != "NA":
                append = True
            if append:
                line = "\t".join(line) + "\n"
                out += line
        f_in.close()

df = pd.read_csv(io.StringIO(out), low_memory=False, sep="\t")

columns = ["target_origin", "query_origin", "overlap_score", 
           "query_synteny_conserved", "target_synteny_conserved", "target_id",
           "query_id"]
ascending = [True, True, False, False, False, True, True]

df = df.astype(str)


df.sort_values(by=columns, ascending=ascending, na_position='last',
               inplace=True)



df.to_csv(final_path + out_file, index=False, sep="\t", na_rep="NA")

out = ""
x = 0
for file in corrected_files:
    if "ref" in file and "filtered" not in file:
        x += 1
        out_file = "all_40X_ref_equivalences.csv"
        j = 0
        f_in = open(out_path + file, "r", encoding="utf-8")
        for line in f_in:
            append = False
            line = line.strip().split("\t")
            j += 1
            if j == 1:
                if x == 1:
                    append = True
            elif line[0] != "NA":
                append = True
            if append:
                line = "\t".join(line) + "\n"
                out += line
        f_in.close()

df = pd.read_csv(io.StringIO(out), low_memory=False, sep="\t")
df.sort_values(by=columns, ascending=ascending, na_position='last',
               inplace=True)
df.to_csv(final_path + out_file, index=False, sep="\t", na_rep="NA")

out = ""
x = 0
for file in corrected_files:
    if "alt" in file and "filtered" in file:
        x += 1
        out_file = "all_40X_alt_equivalences_filtered.csv"
        j = 0
        f_in = open(out_path + file, "r", encoding="utf-8")
        for line in f_in:
            append = False
            line = line.strip().split("\t")
            j += 1
            if j == 1:
                if x == 1:
                    append = True
            elif line[0] != "NA":
                append = True
            if append:
                line = "\t".join(line) + "\n"
                out += line
        f_in.close()

df = pd.read_csv(io.StringIO(out), low_memory=False, sep="\t")
df.sort_values(by=columns, ascending=ascending, na_position='last',
               inplace=True)
df.to_csv(final_path + out_file, index=False, sep="\t", na_rep="NA")

out = ""
x = 0
for file in corrected_files:
    if "alt" in file and "filtered" not in file:
        x += 1
        out_file = "all_40X_alt_equivalences.csv"
        j = 0
        f_in = open(out_path + file, "r", encoding="utf-8")
        for line in f_in:
            append = False
            line = line.strip().split("\t")
            j += 1
            if j == 1:
                if x == 1:
                    append = True
            elif line[0] != "NA":
                append = True
            if append:
                line = "\t".join(line) + "\n"
                out += line
        f_in.close()

df = pd.read_csv(io.StringIO(out), low_memory=False, sep="\t")
df.sort_values(by=columns, ascending=ascending, na_position='last',
               inplace=True)
df.to_csv(final_path + out_file, index=False, sep="\t", na_rep="NA")