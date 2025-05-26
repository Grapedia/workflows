#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Creating all annotation objects for a single species.

@author: David Navarro
"""

# INPUT: change for other species
from config.PN40024 import genome_files, annotation_files, annotation_transfer_files, features_path, pickle_path, species

import time
import os
from modules.tools import pickle_load, pickle_save
from modules.genome import Genome
from modules.annotation import Annotation
from modules.plots import pie_chart

start = time.time()
os.system(f"mkdir -p {features_path}")
os.system(f"mkdir -p {pickle_path}")

genomes = {}

for k, v in genome_files.items():
    genomes[k] = Genome(k, v)

k = "v3_from_12Xv2_raw_on_40X_ref"
location = f"{pickle_path}{k}_annotation.pkl"
if os.path.isfile(location):
    print(f"Loading {k} annotation object\n")
    a = pickle_load(location)
else:
    print(f"{k} annotation does not exist")
    annot_name, genome = k.split("_on_")
    a = Annotation(annot_name, annotation_files[k], genomes[genome])
    a.update_stats(export=True, genome=genomes[a.genome])
    pickle_save(location, a)

location = f"{pickle_path}catalogue_v3_on_40X_ref_annotation.pkl"
if os.path.isfile(location):
    a1 = pickle_load(location)
else:
    root = a.path.split("40X_annotation_transfer")[0]
    a1 = a.copy()
    a1.add_gene_names_and_descriptors(file_path=f"../../genomes_and_annotation/grapevine/PN40024/12Xv2_annotation/catalogue/vvi_lab_catalogue_V4_TM39_unofficial_alternative_stats.csv")
    a1.remove_noncoding_transcripts_from_coding_genes()
    a1.remove_genes_without_names()
    a1.id = "catalogue_v3_on_40X_ref"
    a1.name = "catalogue_v3"
    a1.export_gff()
    a1.update_stats(export=True, genome=genomes[a1.genome])
    pickle_save(location, a1)

location = f"{pickle_path}catalogue_v3_pseudo_on_40X_ref_annotation.pkl"
if os.path.isfile(location):
    a_pseudo = pickle_load(location)
else:
    root = a.path.split("40X_annotation_transfer")[0]
    a_pseudo = a.copy()
    a_pseudo.add_gene_names_and_descriptors(file_path=f"../../genomes_and_annotation/grapevine/PN40024/12Xv2_annotation/catalogue/vvi_lab_catalogue_V4_TM39_unofficial_alternative_pseudo.csv")
    a_pseudo.remove_noncoding_transcripts_from_coding_genes()
    a_pseudo.remove_genes_without_names()
    a_pseudo.id = "catalogue_v3_pseudo_on_40X_ref"
    a_pseudo.name = "catalogue_v3_pseudo"
    a_pseudo.export_gff()
    a_pseudo.update_stats(export=True, genome=genomes[a_pseudo.genome])
    pickle_save(location, a_pseudo)

location = f"{pickle_path}catalogue_v3_expressed_on_40X_ref_annotation.pkl"
if os.path.isfile(location):
    a_expressed = pickle_load(location)
else:
    root = a.path.split("40X_annotation_transfer")[0]
    a_expressed = a.copy()
    a_expressed.add_gene_names_and_descriptors(file_path=f"../../genomes_and_annotation/grapevine/PN40024/12Xv2_annotation/catalogue/vvi_lab_catalogue_V4_TM39_unofficial_alternative_expressed.csv")
    a_expressed.remove_noncoding_transcripts_from_coding_genes()
    a_expressed.remove_genes_without_names()
    a_expressed.id = "catalogue_v3_expressed_on_40X_ref"
    a_expressed.name = "catalogue_v3_expressed"
    a_expressed.export_gff()
    a_expressed.update_stats(export=True, genome=genomes[a_expressed.genome])
    pickle_save(location, a_expressed)

objects = [a, a1, a_expressed, a_pseudo]

for o in objects:
    gene_count = 0
    for genes in o.chrs.values():
        for g in genes.values():
            gene_count += 1
    print(f"{o.id} has {gene_count} genes")

k = "final_reduced_on_40X_ref"

location = f"{pickle_path}{k}_annotation.pkl"
if os.path.isfile(location):
    print(f"Loading {k} annotation object\n")
    b = pickle_load(location)
else:
    print(f"{k} annotation does not exist")
    annot_name, genome = k.split("_on_")
    b = Annotation(annot_name, annotation_files[k], genomes[genome])
    b.update_stats(export=True, genome=genomes[b.genome])
    pickle_save(location, b)

a1.detect_gene_overlaps(b)
a1.export_equivalences(stringent=False, return_df=False, NAs=False, export_csv=True)

d = {}
for i in range(12):
    d[str(i)] = 0

d_labels = {}
for i in range(12):
    d_labels[str(i)] = []

for genes in a1.chrs.values():
    for g in genes.values():
        max_score = 0
        for o in g.overlaps["other"]:
            if o.score > max_score:
                max_score = o.score
        d[str(max_score)] += 1
        d_labels[str(max_score)].append(g.id)

for key in d_labels:
    new_label = ""
    count = 0
    for n, g_id in enumerate(d_labels[key]):
        count += 1
        if count <= 10:
            if n == 0 and count == 1:
                new_label += g_id
            else:
                new_label += f";{g_id}"
        else:
            count = 1
            new_label += f"<br>{g_id}"
    d_labels[key] = new_label

pie_chart(list(d.keys()), list(d.values()), f"{b.path}/out_stats/", f"{a1.id}_{b.id}_found_catalogue_genes", f"Overlap Scores of {a1.id} with {b.id} annotation", hovertext_labels=list(d_labels.values()), colours="extreme")


a_pseudo.detect_gene_overlaps(b)
a_pseudo.export_equivalences(stringent=False, return_df=False, NAs=False, export_csv=True)

d = {}
for i in range(12):
    d[str(i)] = 0

d_labels = {}
for i in range(12):
    d_labels[str(i)] = []

for genes in a_pseudo.chrs.values():
    for g in genes.values():
        max_score = 0
        for o in g.overlaps["other"]:
            if o.score > max_score:
                max_score = o.score
        d[str(max_score)] += 1
        d_labels[str(max_score)].append(g.id)

for key in d_labels:
    new_label = ""
    count = 0
    for n, g_id in enumerate(d_labels[key]):
        count += 1
        if count <= 10:
            if n == 0 and count == 1:
                new_label += g_id
            else:
                new_label += f";{g_id}"
        else:
            count = 1
            new_label += f"<br>{g_id}"
    d_labels[key] = new_label

pie_chart(list(d.keys()), list(d.values()), f"{b.path}/out_stats/", f"{a_pseudo.id}_{b.id}_found_catalogue_genes", f"Overlap Scores of {a_pseudo.id} with {b.id} annotation", hovertext_labels=list(d_labels.values()), colours="extreme")

a_expressed.detect_gene_overlaps(b)
a_expressed.export_equivalences(stringent=False, return_df=False, NAs=False, export_csv=True)

d = {}
for i in range(12):
    d[str(i)] = 0

d_labels = {}
for i in range(12):
    d_labels[str(i)] = []

for genes in a_expressed.chrs.values():
    for g in genes.values():
        max_score = 0
        for o in g.overlaps["other"]:
            if o.score > max_score:
                max_score = o.score
        d[str(max_score)] += 1
        d_labels[str(max_score)].append(g.id)

for key in d_labels:
    new_label = ""
    count = 0
    for n, g_id in enumerate(d_labels[key]):
        count += 1
        if count <= 10:
            if n == 0 and count == 1:
                new_label += g_id
            else:
                new_label += f";{g_id}"
        else:
            count = 1
            new_label += f"<br>{g_id}"
    d_labels[key] = new_label

pie_chart(list(d.keys()), list(d.values()), f"{b.path}/out_stats/", f"{a_expressed.id}_{b.id}_found_catalogue_genes", f"Overlap Scores of {a_expressed.id} with {b.id} annotation", hovertext_labels=list(d_labels.values()), colours="extreme")

a1.clear_overlaps()
a_pseudo.clear_overlaps()
a_expressed.clear_overlaps()

k = "v4.1_on_40X_ref"

location = f"{pickle_path}{k}_annotation.pkl"
if os.path.isfile(location):
    print(f"Loading {k} annotation object\n")
    b = pickle_load(location)
else:
    print(f"{k} annotation does not exist")
    annot_name, genome = k.split("_on_")
    b = Annotation(annot_name, annotation_files[k], genomes[genome])
    b.update_stats(export=True, genome=genomes[b.genome])
    pickle_save(location, b)

a1.detect_gene_overlaps(b)
a1.export_equivalences(stringent=False, return_df=False, NAs=False, export_csv=True)

d = {}
for i in range(12):
    d[str(i)] = 0

d_labels = {}
for i in range(12):
    d_labels[str(i)] = []

for genes in a1.chrs.values():
    for g in genes.values():
        max_score = 0
        for o in g.overlaps["other"]:
            if o.score > max_score:
                max_score = o.score
        d[str(max_score)] += 1
        d_labels[str(max_score)].append(g.id)

for key in d_labels:
    new_label = ""
    count = 0
    for n, g_id in enumerate(d_labels[key]):
        count += 1
        if count <= 10:
            if n == 0 and count == 1:
                new_label += g_id
            else:
                new_label += f";{g_id}"
        else:
            count = 1
            new_label += f"<br>{g_id}"
    d_labels[key] = new_label

pie_chart(list(d.keys()), list(d.values()), f"{b.path}/out_stats/", f"{a1.id}_{b.id}_found_catalogue_genes", f"Overlap Scores of {a1.id} with {b.id} annotation", hovertext_labels=list(d_labels.values()), colours="extreme")


a_pseudo.detect_gene_overlaps(b)
a_pseudo.export_equivalences(stringent=False, return_df=False, NAs=False, export_csv=True)

d = {}
for i in range(12):
    d[str(i)] = 0

d_labels = {}
for i in range(12):
    d_labels[str(i)] = []

for genes in a_pseudo.chrs.values():
    for g in genes.values():
        max_score = 0
        for o in g.overlaps["other"]:
            if o.score > max_score:
                max_score = o.score
        d[str(max_score)] += 1
        d_labels[str(max_score)].append(g.id)

for key in d_labels:
    new_label = ""
    count = 0
    for n, g_id in enumerate(d_labels[key]):
        count += 1
        if count <= 10:
            if n == 0 and count == 1:
                new_label += g_id
            else:
                new_label += f";{g_id}"
        else:
            count = 1
            new_label += f"<br>{g_id}"
    d_labels[key] = new_label

pie_chart(list(d.keys()), list(d.values()), f"{b.path}/out_stats/", f"{a_pseudo.id}_{b.id}_found_catalogue_genes", f"Overlap Scores of {a_pseudo.id} with {b.id} annotation", hovertext_labels=list(d_labels.values()), colours="extreme")

a_expressed.detect_gene_overlaps(b)
a_expressed.export_equivalences(stringent=False, return_df=False, NAs=False, export_csv=True)

d = {}
for i in range(12):
    d[str(i)] = 0

d_labels = {}
for i in range(12):
    d_labels[str(i)] = []

for genes in a_expressed.chrs.values():
    for g in genes.values():
        max_score = 0
        for o in g.overlaps["other"]:
            if o.score > max_score:
                max_score = o.score
        d[str(max_score)] += 1
        d_labels[str(max_score)].append(g.id)

for key in d_labels:
    new_label = ""
    count = 0
    for n, g_id in enumerate(d_labels[key]):
        count += 1
        if count <= 10:
            if n == 0 and count == 1:
                new_label += g_id
            else:
                new_label += f";{g_id}"
        else:
            count = 1
            new_label += f"<br>{g_id}"
    d_labels[key] = new_label

pie_chart(list(d.keys()), list(d.values()), f"{b.path}/out_stats/", f"{a_expressed.id}_{b.id}_found_catalogue_genes", f"Overlap Scores of {a_expressed.id} with {b.id} annotation", hovertext_labels=list(d_labels.values()), colours="extreme")

a1.clear_overlaps()
a_pseudo.clear_overlaps()
a_expressed.clear_overlaps()

k = "EVM_PASA_on_40X_ref"

location = f"{pickle_path}{k}_annotation.pkl"
if os.path.isfile(location):
    print(f"Loading {k} annotation object\n")
    b = pickle_load(location)
else:
    print(f"{k} annotation does not exist")
    annot_name, genome = k.split("_on_")
    b = Annotation(annot_name, annotation_files[k], genomes[genome])
    b.update_stats(export=True, genome=genomes[b.genome])
    pickle_save(location, b)

a1.detect_gene_overlaps(b)
a1.export_equivalences(stringent=False, return_df=False, NAs=False, export_csv=True)

d = {}
for i in range(12):
    d[str(i)] = 0

d_labels = {}
for i in range(12):
    d_labels[str(i)] = []

for genes in a1.chrs.values():
    for g in genes.values():
        max_score = 0
        for o in g.overlaps["other"]:
            if o.score > max_score:
                max_score = o.score
        d[str(max_score)] += 1
        d_labels[str(max_score)].append(g.id)

for key in d_labels:
    new_label = ""
    count = 0
    for n, g_id in enumerate(d_labels[key]):
        count += 1
        if count <= 10:
            if n == 0 and count == 1:
                new_label += g_id
            else:
                new_label += f";{g_id}"
        else:
            count = 1
            new_label += f"<br>{g_id}"
    d_labels[key] = new_label

pie_chart(list(d.keys()), list(d.values()), f"{b.path}/out_stats/", f"{a1.id}_{b.id}_found_catalogue_genes", f"Overlap Scores of {a1.id} with {b.id} annotation", hovertext_labels=list(d_labels.values()), colours="extreme")


a_pseudo.detect_gene_overlaps(b)
a_pseudo.export_equivalences(stringent=False, return_df=False, NAs=False, export_csv=True)

d = {}
for i in range(12):
    d[str(i)] = 0

d_labels = {}
for i in range(12):
    d_labels[str(i)] = []

for genes in a_pseudo.chrs.values():
    for g in genes.values():
        max_score = 0
        for o in g.overlaps["other"]:
            if o.score > max_score:
                max_score = o.score
        d[str(max_score)] += 1
        d_labels[str(max_score)].append(g.id)

for key in d_labels:
    new_label = ""
    count = 0
    for n, g_id in enumerate(d_labels[key]):
        count += 1
        if count <= 10:
            if n == 0 and count == 1:
                new_label += g_id
            else:
                new_label += f";{g_id}"
        else:
            count = 1
            new_label += f"<br>{g_id}"
    d_labels[key] = new_label

pie_chart(list(d.keys()), list(d.values()), f"{b.path}/out_stats/", f"{a_pseudo.id}_{b.id}_found_catalogue_genes", f"Overlap Scores of {a_pseudo.id} with {b.id} annotation", hovertext_labels=list(d_labels.values()), colours="extreme")

a_expressed.detect_gene_overlaps(b)
a_expressed.export_equivalences(stringent=False, return_df=False, NAs=False, export_csv=True)

d = {}
for i in range(12):
    d[str(i)] = 0

d_labels = {}
for i in range(12):
    d_labels[str(i)] = []

for genes in a_expressed.chrs.values():
    for g in genes.values():
        max_score = 0
        for o in g.overlaps["other"]:
            if o.score > max_score:
                max_score = o.score
        d[str(max_score)] += 1
        d_labels[str(max_score)].append(g.id)

for key in d_labels:
    new_label = ""
    count = 0
    for n, g_id in enumerate(d_labels[key]):
        count += 1
        if count <= 10:
            if n == 0 and count == 1:
                new_label += g_id
            else:
                new_label += f";{g_id}"
        else:
            count = 1
            new_label += f"<br>{g_id}"
    d_labels[key] = new_label

pie_chart(list(d.keys()), list(d.values()), f"{b.path}/out_stats/", f"{a_expressed.id}_{b.id}_found_catalogue_genes", f"Overlap Scores of {a_expressed.id} with {b.id} annotation", hovertext_labels=list(d_labels.values()), colours="extreme")

now = time.time()
lapse = now - start
print(f"\nCreating and updating {species} annotation objects as well as "
      f"exporting all features and gene lengths took {round((lapse/60)/60, 1)} "
       "hours\n")


