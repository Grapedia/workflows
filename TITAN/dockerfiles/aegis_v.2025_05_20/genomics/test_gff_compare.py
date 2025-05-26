#!/usr/bin/env python3
# -*- coding: utf-8 -*-
root = "../../genomes_and_annotation/grapevine/PN40024/40X_annotation/"
f1 = f"{root}final_annotation_reduced_redundancy_on_40X_test9_output.gff3"
f2 = f"{root}out_gffs/test9_export_b4_renaming_on_40X_ref_output.gff3"
f3 = f"{root}out_gffs/test9_export_on_40X_ref_output.gff3"

common = ["gene", "mRNA", "exon", "ncRNA", "CDS"]

f1_d = {}
for c in common:
    f1_d[c] = {}
f1_d["other"] = {}

f_in = open(f1, encoding="utf-8")
for line in f_in:
    if line.startswith("#"):
        continue
    line = line.strip().split("\t")
    ft = line[2]
    attributes = line[-1]
    start = line[3]
    end = line[4]
    ch = line[0]
    id = ""
    parents = []
    for a in attributes.split(";"):
        if a.split("=")[0] == "ID":
            id = a.split("=")[1].strip()
        elif a.split("=")[0] == "Parent":
            parents = a.split("=")[1].strip().split(",")
    id = f"{id}_{start}_{end}"
    if ft in common:
        if id not in f1_d[ft]:
            f1_d[ft][id] = (ch, start, end, ch, parents)
        else:
            print(f"f1 id problem: {id} {ft}")
    else:
        if id not in f1_d["other"]:
            f1_d["other"][id] = (ch, start, end, ch, parents)
        else:
            print(f"f1 id problem: {id} {ft}")
f_in.close()

f2_d = {}
for c in common:
    f2_d[c] = {}
f2_d["other"] = {}

f_in = open(f2, encoding="utf-8")
for line in f_in:
    if line.startswith("#"):
        continue
    line = line.strip().split("\t")
    ft = line[2]
    attributes = line[-1]
    start = line[3]
    end = line[4]
    ch = line[0]
    id = ""
    parents = []
    for a in attributes.split(";"):
        if a.split("=")[0] == "ID":
            id = a.split("=")[1].strip()
        elif a.split("=")[0] == "Parent":
            parents = a.split("=")[1].strip().split(",")
    id = f"{id}_{start}_{end}"
    if ft in common:
        if id not in f2_d[ft]:
            f2_d[ft][id] = (ch, start, end, ch, parents)
        else:
            print(f"f2 id problem: {id} {ft}")
    else:
        if id not in f2_d["other"]:
            f2_d["other"][id] = (ch, start, end, ch, parents)
        else:
            print(f"f2 id problem: {id} {ft}")
f_in.close()

f3_d = {}
for c in common:
    f3_d[c] = {}
f3_d["other"] = {}

f_in = open(f3, encoding="utf-8")
for line in f_in:
    if line.startswith("#"):
        continue
    line = line.strip().split("\t")
    ft = line[2]
    attributes = line[-1]
    start = line[3]
    end = line[4]
    ch = line[0]
    id = ""
    parents = []
    for a in attributes.split(";"):
        if a.split("=")[0] == "ID":
            id = a.split("=")[1].strip()
        elif a.split("=")[0] == "Parent":
            parents = a.split("=")[1].strip().split(",")
    id = f"{id}_{start}_{end}"
    if ft in common:
        if id not in f3_d[ft]:
            f3_d[ft][id] = (ch, start, end, ch, parents)
        else:
            print(f"f3 id problem: {id} {ft}")
    else:
        if id not in f3_d["other"]:
            f3_d["other"][id] = (ch, start, end, ch, parents)
        else:
            print(f"f3 id problem: {id} {ft}")
f_in.close()

tags = ["input", "output", "output_renamed"]
ds = [f1_d, f2_d, f3_d]

for n, d in enumerate(ds):
    for ft, fts in d.items():
        print(f"There are {len(fts.keys())} {ft} features in {tags[n]}")

for id in f1_d["exon"]:
    if id not in f2_d["exon"]:
        print(id)

for id in f1_d["CDS"]:
    if id not in f2_d["CDS"]:
        print(id)

for id in f2_d["CDS"]:
    if id not in f1_d["CDS"]:
        print(id)