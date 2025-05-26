#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Comparing two gffs in terms of gene, mRNA, exon, lncRNA coordinates to see if they are exactly the same,
ignores UTRs.

@author: David Navarro
"""

from deepdiff import DeepDiff

f1 = "../../genomes_and_annotation/grapevine/PN40024/test_gffs/5.1_on_T2T_ref_new.gff3"
f2 = "../../genomes_and_annotation/grapevine/PN40024/test_gffs/5.1_on_T2T_ref_new2.gff3"


files = [f1, f2]

dictionaries = [{}, {}]


for x, file in enumerate(files):
    print(file)
    f_in = open(file)
    for line in f_in:
        if line.startswith("#"):
            continue
        line = line.strip().split("\t")
        attributes = line[-1].split(";")
        ID = None
        start = line[3]
        end = line[4]
        for a in attributes:
            if a.startswith("ID="):
                ID = a.split("=")[1].strip()
        feature = line[2]
        if "UTR" in feature:
            continue

        if feature not in dictionaries[x]:
            dictionaries[x][feature] = {}
        
        if ID not in dictionaries[x][feature]:
            dictionaries[x][feature][ID] = (start, end)

        elif feature == "CDS":
            new = [dictionaries[x][feature][ID][0], dictionaries[x][feature][ID][1]]
            if start < dictionaries[x][feature][ID][0]:
                new[0] = start
            if end > dictionaries[x][feature][ID][1]:
                new[1] = end
            dictionaries[x][feature][ID] = tuple(new)
        else:
            print(f"Error: there should not be more than one {ID}")

    f_in.close()

diff = DeepDiff(dictionaries[0], dictionaries[1])


print(diff)

        


if dictionaries[0] != dictionaries[1]:
    print("Gffs are different in terms of gene, mRNA, lncRNA, CDS or exon coordinates and ids")
else:
    print("Gffs are identical in terms of gene, mRNA, lncRNA, CDS or exon coordinates and ids")