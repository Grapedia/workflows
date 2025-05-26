#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
2023/05/12

Program that transforms a gff into a gff compatible with the exon-based
featureCounts summarising parameter required by the SRA_transcriptomics
repository.

This program assumes:
- that the highest level feature in the gff is a gene
- that all the features after a particular gene belong to that gene and not a
different one, i.e. that the gff features are in correct order after gt sort

@author: David Navarro
"""

import os

input_gffs_folder = "input_gffs/"
tidy_gffs_folder = "output_tidy_gffs/"
tidy_fc_gffs_folder = "output_tidy_fc_gffs/"
tidy_fc_simple_gffs_folder = "output_tidy_fc_simple_gffs/"

os.system(f"mkdir -p {input_gffs_folder}")
os.system(f"touch ./{input_gffs_folder}.gitkeep")
os.system(f"mkdir -p {tidy_gffs_folder}")
os.system(f"mkdir -p {tidy_fc_gffs_folder}")
os.system(f"mkdir -p {tidy_fc_simple_gffs_folder}")

input_gffs = os.listdir(input_gffs_folder)
# output extension
extension=".gff3"

for file in input_gffs:
    if not ".gff" in file:
        continue
    file = file.split(".")[0:-1]
    file = ".".join(file)
    tag = file + "_tidy"
    gt_command = f"gt gff3 -sort -tidy -retainids -o {tidy_gffs_folder}{tag}{extension} {input_gffs_folder}{file}{extension}"
    
    if os.path.isfile(tidy_gffs_folder + tag + extension):
        continue
    print(f"Running genome tools gff3 -sort, -tidy and -retainids for {file}")
    os.system(gt_command)


tidy_gffs = os.listdir(tidy_gffs_folder)

# features that should not carry a featurecounts_id since they are not 
# associated to genes
gene_independent_features = ["region", "sequence_feature", "cDNA_match", 
                             "match", "chromosome", "transcript_region",
                             "transposon_fragment"]

# assumption: features are ordered after the gt command above
# second assumption: tricky features that are not to be assigned to genes
# are in the variable above, check for new ones in new gff files
for file in tidy_gffs:
    if not ".gff" in file:
        continue
    file_tag = file.split(".")[0:-1]
    file_tag = ".".join(file_tag)
    file_out = tidy_fc_gffs_folder + file_tag + "_fc" + extension
    out = ""
    f_in = open(tidy_gffs_folder + file, "r", encoding="utf-8")
    for line in f_in:
        if not line or line.startswith("#") or line.startswith("@"):
            out += line
        else:
            line = line.strip()
            temp = line.split("\t")
            if temp[2] == "gene" or temp[2] == "pseudogene" or temp[2] == "transposable_element_gene":
                try:
                    gene_id = temp[8].split(";")[0].split("=")[1]
                except:
                    print(f"file '{file}' with line '{temp}' split failed")
            
            if temp[2] not in gene_independent_features:
                try:
                    line += f";featurecounts_id={gene_id}"
                except:
                    print(f"file '{file}' with line '{temp}' appending failed")
            line += "\n"
            out += line
    f_in.close()
    if os.path.isfile(file_out):
        continue
    print(f"Formatting {file} following featureCounts "
          "exon summary requirements")
    f_out = open(file_out, "w", encoding="utf-8")
    f_out.write(out)
    f_out.close()

fc_gffs = os.listdir(tidy_fc_gffs_folder)
description_fields_to_remove = ["curator_summary", "nochangenat-description",
                                "computational_description", "description",
                                "Note", "Name", "Type", "Prot_size", "Genoscope_ID", "CRIBI_V1_ID"]
# if not empty it trumps above
#description_fields_to_keep = []
description_fields_to_keep = ["ID", "Parent"]

for file in fc_gffs:
    if not ".gff" in file:
        continue
    file_tag = file.split(".")[0:-1]
    file_tag = ".".join(file_tag)
    file_out = tidy_fc_simple_gffs_folder + file_tag + "_simple" + extension
    out = ""
    f_in = open(tidy_fc_gffs_folder + file, "r", encoding="utf-8")
    for line in f_in:
        if not line or line.startswith("#") or line.startswith("@"):
            out += line
        else:
            line = line.strip()
            temp = line.split("\t")
            temp[-1] = temp[-1].split(";")
            new_atts = []
            for a in temp[-1]:
                term = a.split("=")[0]
                if description_fields_to_keep != []:
                    if term in description_fields_to_keep:
                        new_atts.append(a)
                else:
                    if term not in description_fields_to_remove:
                        new_atts.append(a)
            temp[-1] = ";".join(new_atts)
            line = "\t".join(temp)
            line += "\n"
            out += line
    f_in.close()
    if os.path.isfile(file_out):
        continue
    print(f"Formatting {file} to remove unwanted description attributes")
    f_out = open(file_out, "w", encoding="utf-8")
    f_out.write(out)
    f_out.close()