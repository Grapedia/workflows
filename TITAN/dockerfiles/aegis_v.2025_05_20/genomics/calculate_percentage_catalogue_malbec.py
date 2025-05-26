#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Creating all annotation objects for a single species in malbec.

"""

import time
import os
from modules.tools import pickle_load, pickle_save
from modules.geneclasses import Genome, Annotation
from modules.plots import pie_chart

start = time.time()


g = Genome('Malbec_Magde','/media/inigo/F4DE8D50DE8D0BD43/inigo/genomes_and_annotation/grapevine/malbec/magdeleine_assembly/Vitis_malbecmagdeleine.Malbec_Magdeleine_v1.0.dna.toplevel.fa')

pickle_liftoff = "/media/inigo/F4DE8D50DE8D0BD43/inigo/annotation/Malbec/genomes/liftoff_catalogue_v3_on_malbec_magde.pkl"

if os.path.isfile(pickle_liftoff):
    print(f"Loading v3 annotation object\n")
    start = time.time()
    liftoff = pickle_load(pickle_liftoff)
    now = time.time()
    lapse = now - start
    print(f"\nCreating v3 annotation object took {round(lapse, 1)} "
              "seconds\n")
else:
    print(f"v3 annotation does not exist\n")
    liftoff = Annotation('v3_malbec_liftoff', '/media/inigo/F4DE8D50DE8D0BD43/inigo/annotation/Malbec/Magdeleine/liftoff/v3/Malbec_magde_v3.gff3', g)
    pickle_save(pickle_liftoff, liftoff)

pickle_malbec = "/media/inigo/F4DE8D50DE8D0BD43/inigo/annotation/Malbec/genomes/malbec_magde_catalogue_v3_on_malbec.pkl"

if os.path.isfile(pickle_malbec):
    print(f"Loading malbec annotation object\n")
    start = time.time()
    malbec = pickle_load(pickle_malbec)
    now = time.time()
    lapse = now - start
    print(f"\nCreating malbec_magde annotation object took {round(lapse, 1)} "
              "seconds\n")
else:
    print(f"Malbec_magde annotation does not exist")
    malbec = Annotation('MalbecMagde_annotation', '/media/inigo/F4DE8D50DE8D0BD43/inigo/annotation/Malbec/genomes/vitis_malbecmagdeleine.gff', g)
    pickle_save(pickle_malbec, malbec)

malbec.path = '/media/inigo/F4DE8D50DE8D0BD43/inigo/annotation/Malbec/genomes/'
liftoff.path = '/media/inigo/F4DE8D50DE8D0BD43/inigo/annotation/Malbec/Prunelard/liftoff/v3/'


liftoff_all = liftoff.copy()
liftoff_pseudo = liftoff.copy()
liftoff_expressed = liftoff.copy()

del liftoff


liftoff_all.add_gene_names_and_descriptors(file_path=f"../../genomes_and_annotation/grapevine/PN40024/12Xv2_annotation/catalogue/vvi_lab_catalogue_V4_TM39_unofficial_alternative_stats.csv")

liftoff_all.remove_genes_without_names()

liftoff_pseudo.add_gene_names_and_descriptors(file_path=f"../../genomes_and_annotation/grapevine/PN40024/12Xv2_annotation/catalogue/vvi_lab_catalogue_V4_TM39_unofficial_alternative_pseudo.csv")

liftoff_pseudo.remove_genes_without_names()

liftoff_expressed.add_gene_names_and_descriptors(file_path=f"../../genomes_and_annotation/grapevine/PN40024/12Xv2_annotation/catalogue/vvi_lab_catalogue_V4_TM39_unofficial_alternative_expressed.csv")

liftoff_expressed.remove_genes_without_names()


liftoff_all.id = "catalogue_v3_on_malbec_magde"
liftoff_all.name = "catalogue_v3"

liftoff_pseudo.id = "catalogue_v3_on_malbec_magde_pseudo"
liftoff_pseudo.name = "catalogue_v3_pseudo"

liftoff_expressed.id = "catalogue_v3_on_malbec_magde_expressed"
liftoff_expressed.name = "catalogue_v3_expressed"


liftoff_all.detect_gene_overlaps(malbec)
liftoff_all.export_equivalences(stringent=False, return_df=False, NAs=False, export_csv=True)


def create_pie_chart_dictionaries(a):
    d = {}
    d["missing"] = 0
    d["present"] = 0

    d_labels = {}
    d_labels["missing"] = []
    d_labels["present"] = []

    for genes in a.chrs.values():
        for g in genes.values():
            max_score = 0
            for o in g.overlaps["other"]:
                if o.score > max_score:
                    max_score = o.score
            if max_score > 5:
                d["present"] += 1
                d_labels["present"].append(g.id)
            else:
                d["missing"] += 1
                d_labels["missing"].append(g.id)

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

    return d, d_labels



def create_pie_chart_full_dictionaries(a):
    d = {}
    for i in range(12):
        d[str(i)] = 0

    d_labels = {}
    for i in range(12):
        d_labels[str(i)] = []

    for genes in a.chrs.values():
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

    return d, d_labels


#Liftoff_all

d, d_labels = create_pie_chart_dictionaries(liftoff_all)

pie_chart(list(d.keys()), list(d.values()), f"{malbec.path}/out_stats/", f"{liftoff_all.id}_{malbec.id}_found_catalogue_genes_simplified", f"Simplified overlap scores of {liftoff_all.id} with {malbec.id} annotation", hovertext_labels=list(d_labels.values()), colours="extreme")


d, d_labels = create_pie_chart_full_dictionaries(liftoff_all)

pie_chart(list(d.keys()), list(d.values()), f"{malbec.path}/out_stats/", f"{liftoff_all.id}_{malbec.id}_found_catalogue_genes", f"Overlap scores of {liftoff_all.id} with {malbec.id} annotation", hovertext_labels=list(d_labels.values()), colours="extreme")

#Liftoff_pseudo

liftoff_pseudo.detect_gene_overlaps(malbec)
liftoff_pseudo.export_equivalences(stringent=False, return_df=False, NAs=False, export_csv=True)

d, d_labels = create_pie_chart_dictionaries(liftoff_pseudo)

pie_chart(list(d.keys()), list(d.values()), f"{malbec.path}/out_stats/", f"{liftoff_pseudo.id}_{malbec.id}_found_catalogue_genes_simplified", f"Simplified overlap Scores of {liftoff_pseudo.id} with {malbec.id} annotation", hovertext_labels=list(d_labels.values()), colours="extreme")

d, d_labels = create_pie_chart_full_dictionaries(liftoff_pseudo)

pie_chart(list(d.keys()), list(d.values()), f"{malbec.path}/out_stats/", f"{liftoff_pseudo.id}_{malbec.id}_found_catalogue_genes", f"Overlap scores of {liftoff_pseudo.id} with {malbec.id} annotation", hovertext_labels=list(d_labels.values()), colours="extreme")

#Liftoff_expressed

liftoff_expressed.detect_gene_overlaps(malbec)
liftoff_expressed.export_equivalences(stringent=False, return_df=False, NAs=False, export_csv=True)

d, d_labels = create_pie_chart_dictionaries(liftoff_expressed)


pie_chart(list(d.keys()), list(d.values()), f"{malbec.path}/out_stats/", f"{liftoff_expressed.id}_{malbec.id}_found_catalogue_genes_simplified", f"Simplified overlap Scores of {liftoff_expressed.id} with {malbec.id} annotation", hovertext_labels=list(d_labels.values()), colours="extreme")

d, d_labels = create_pie_chart_full_dictionaries(liftoff_expressed)

pie_chart(list(d.keys()), list(d.values()), f"{malbec.path}/out_stats/", f"{liftoff_expressed.id}_{malbec.id}_found_catalogue_genes", f"Overlap scores of {liftoff_expressed.id} with {malbec.id} annotation", hovertext_labels=list(d_labels.values()), colours="extreme")


























