#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Creating all annotation objects for a single species.

@author: David Navarro
"""

# INPUT: change for other species
from config.PN40024_new_T2T import genome_files, pre_annotation_files, features_path, pickle_path, species

import time
import os
from modules.tools import pickle_load, pickle_save
from modules.genome import Genome
from modules.annotation import Annotation

start = time.time()
os.system(f"mkdir -p {features_path}")
os.system(f"mkdir -p {pickle_path}")

genomes = {}

for k, v in genome_files.items():
    genomes[k] = Genome(k, v)
    genomes[k].export_feature_sizes()

annotations = {}

for k, v in pre_annotation_files.items():
    location = f"{pickle_path}{k}_annotation.pkl"
    if os.path.isfile(location):
        print(f"Loading {k} annotation object\n")
        annotations[k] = pickle_load(location)
    else:
        annot_name, genome = k.split("_on_")
        annotations[k] = Annotation(annot_name, v, genomes[genome], rework_CDSs=False)
        pickle_save(location, annotations[k])


# add new non-overlapping genes (including fully intron nested), i.e. fully non-coding genes!

a = annotations["v5.1_ncRNAs_removed_on_T2T_ref"].copy()

num_transcripts = 0
num_genes = 0
for genes in a.chrs.values():
    for g in genes.values():
        num_genes += 1
        for t in g.transcripts.values():
            num_transcripts += 1

print(f"There are {num_genes} genes and {num_transcripts} transcripts in the original annotation.")

num_transcripts = 0
num_ncRNAs = 0
num_genes = 0
for genes in annotations["just_lncRNAs_on_T2T_ref"].chrs.values():
    num_genes += len(genes)
    for g in genes.values():
        num_transcripts += 1
        for t in g.transcripts.values():
            num_ncRNAs += 1
            t.feature = "lncRNA"
            if t.attributes != "":
                t.attributes += f";nc_stringtie={g.id}"
            else:
                t.attributes = f"nc_stringtie={g.id}"

print(f"There are {num_genes} ncRNA genes and {num_ncRNAs} transcripts out of a total of {num_transcripts} detected by the lncRNA pipeline.")

a.merge(annotations["just_lncRNAs_on_T2T_ref"], ignore_overlaps=False, exon_overlap_threshold=0)

num_ncRNAs = 0
num_genes = 0
num_transcripts = 0
for genes in a.chrs.values():
    for g in genes.values():
        num_genes += 1
        for t in g.transcripts.values():
            num_transcripts += 1
            if t.feature == "lncRNA":
                num_ncRNAs += 1

print(f"There are {num_ncRNAs} ncRNA transcripts out of a total {num_transcripts} and {num_genes} genes after the first merge.")

# adding lncRNAs to coding genes as transcript variants

a.detect_gene_overlaps(annotations["just_lncRNAs_on_T2T_ref"])

already_added_genes = {}

for chrom, genes in annotations["just_lncRNAs_on_T2T_ref"].chrs.items():
    for g in genes.values():
        for o in g.overlaps["other"]:
            if o.exon_query_percent == 100:
                lncRNA = False
                for t in a.chrs[chrom][o.id].transcripts.values():
                    if t.feature == "lncRNA":
                        lncRNA = True
                if lncRNA:
                    already_added_genes[g.id] = chrom

for g_id, chrom in already_added_genes.items():
    del annotations["just_lncRNAs_on_T2T_ref"].chrs[chrom][g_id]

annotations["just_lncRNAs_on_T2T_ref"].update()

num_ncRNAs = 0
num_genes = 0
num_transcripts = 0
for genes in annotations["just_lncRNAs_on_T2T_ref"].chrs.values():
    for g in genes.values():
        num_genes += 1
        for t in g.transcripts.values():
            num_transcripts += 1
            if t.feature == "lncRNA":
                num_ncRNAs += 1

print(f"There are {num_ncRNAs} ncRNA transcripts out of a total {num_transcripts} and {num_genes} genes remaining in the potential ncRNAs annotation.")

a.detect_gene_overlaps(annotations["just_lncRNAs_on_T2T_ref"])

for chrom, genes in annotations["just_lncRNAs_on_T2T_ref"].chrs.items():
    for g in genes.values():
        max_overlap_value = 0
        max_overlap_gene = None
        for o in g.overlaps["other"]:
            if o.exon_query_percent > 0:
                if o.exon_query_percent > max_overlap_value:
                    max_overlap_value = o.exon_query_percent
                    max_overlap_gene = o.id
        if max_overlap_gene != None:
            counter = len(a.chrs[chrom][max_overlap_gene].transcripts)
            for t in g.transcripts.values():
                counter += 1
                temp = t.copy()
                temp.id = f"{max_overlap_gene}_t_extra_{counter}"

                if a.chrs[chrom][max_overlap_gene].strand == temp.strand:
                    if t.attributes != "":
                        t.attributes += f";antisense=False"
                    else:
                        t.attributes = f"antisense=False"
                else:
                    if t.attributes != "":
                        t.attributes += f";antisense=True"
                    else:
                        t.attributes = f"antisense=True"

                a.chrs[chrom][max_overlap_gene].transcripts[f"{max_overlap_gene}_t_extra_{counter}"] = temp
                a.all_transcript_ids[f"{max_overlap_gene}_t_extra_{counter}"] = chrom

num_ncRNAs = 0
num_genes = 0
num_transcripts = 0
for genes in a.chrs.values():
    for g in genes.values():
        num_genes += 1
        for t in g.transcripts.values():
            num_transcripts += 1
            if t.feature == "lncRNA":
                num_ncRNAs += 1

print(f"There are {num_ncRNAs} ncRNA transcripts out of a total {num_transcripts} and {num_genes} genes in the final pre-release annotation.")

id_prefix = "Vitvi05_01"
source = "Titan"
name = "5.1"
id = "5.1_on_T2T_ref"
a.release(custom_path=a.path, name=name, id=id, source_name=source, id_prefix=id_prefix, spacer=10)


a.name = "5.1_coding_only"
a.id = "5.1_coding_only_on_T2T_ref"


a.export_gff(custom_path=a.path, UTRs=True, exclude_non_coding=True)

num_ncRNAs = 0
num_genes = 0
num_transcripts = 0
for genes in a.chrs.values():
    for g in genes.values():
        num_genes += 1
        for t in g.transcripts.values():
            num_transcripts += 1
            if t.feature == "lncRNA":
                num_ncRNAs += 1

print(f"There are {num_ncRNAs} ncRNA transcripts out of a total {num_transcripts} and {num_genes} genes in the final release annotation.")

now = time.time()
lapse = now - start
print(f"\nCreating and updating {species} annotation objects as well as exporting all features and gene lengths took {round((lapse/60)/60, 1)} hours\n")