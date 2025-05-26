#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Calculating percentage presence of catalogue genes (grapevine only)

@author: David Navarro
"""

# INPUT: change for other species
from config.grapevine_malbec import genome_files, annotation_files, features_path, pickle_path, species

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

def create_pie_chart_dictionaries(annotation, detailed:bool=False):
    d = {}
    d_labels = {}

    if not detailed:
        d["missing"] = 0
        d["present"] = 0
        
        d_labels["missing"] = []
        d_labels["present"] = []

        for genes in annotation.chrs.values():
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
    else:
        for i in range(12):
            d[str(i)] = 0

        for i in range(12):
            d_labels[str(i)] = []

        for genes in annotation.chrs.values():
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

def plot_catalogue_stats(query_annotation_name, catalogue_annotation_name, catalogue_name, catalogue_path, catalogue_expressed_path=None, catalogue_pseudo_path=None, chosen_chromosomes:list=None):

    location = f"{pickle_path}{query_annotation_name}_annotation.pkl"
    if os.path.isfile(location):
        print(f"Loading {query_annotation_name} annotation object\n")
        query_a = pickle_load(location)
    else:
        print(f"{query_annotation_name} annotation does not exist")
        annot_name, genome = query_annotation_name.split("_on_")
        query_a = Annotation(annot_name, annotation_files[query_annotation_name], genomes[genome], chosen_chromosomes=chosen_chromosomes)
        query_a.update_stats(export=True, genome=genomes[query_a.genome])
        pickle_save(location, query_a)

    location = f"{pickle_path}{catalogue_annotation_name}_annotation.pkl"
    if os.path.isfile(location):
        print(f"Loading {catalogue_annotation_name} annotation object\n")
        catalogue_pre = pickle_load(location)
    else:
        print(f"{catalogue_annotation_name} annotation does not exist")
        annot_name, genome = catalogue_annotation_name.split("_on_")
        catalogue_pre = Annotation(annot_name, annotation_files[catalogue_annotation_name], genomes[genome],chosen_chromosomes=chosen_chromosomes)
        catalogue_pre.update_stats(export=True, genome=genomes[catalogue_pre.genome])
        pickle_save(location, catalogue_pre)

    location = f"{pickle_path}full_{catalogue_name}_annotation.pkl"
    if os.path.isfile(location):
        catalogue_a = pickle_load(location)
    else:
        catalogue_a = catalogue_pre.copy()
        catalogue_a.add_gene_names_and_descriptors(file_path=catalogue_path)
        catalogue_a.remove_noncoding_transcripts_from_coding_genes()
        catalogue_a.remove_genes_without_names()
        catalogue_a.id = f"full_{catalogue_name}"
        catalogue_a.name = catalogue_a.id.split("_on_")[0]
        catalogue_a.export_gff()
        catalogue_a.update_stats(export=True, genome=genomes[catalogue_a.genome])
        pickle_save(location, catalogue_a)

    query_a.clear_overlaps()
    catalogue_a.clear_overlaps()
    catalogue_a.detect_gene_overlaps(query_a)
    catalogue_a.export_equivalences(stringent=False, return_df=False, NAs=False, export_csv=True)

    d, d_labels = create_pie_chart_dictionaries(catalogue_a)

    os.system(f"mkdir -p {query_a.path}/out_stats/")

    pie_chart(list(d.keys()), list(d.values()), f"{query_a.path}/out_stats/", f"{catalogue_a.id}_{query_a.id}_found_catalogue_genes_simplified", f"Simplified overlap scores of {catalogue_a.id} with {query_a.id} annotation", hovertext_labels=list(d_labels.values()), colours="extreme")

    d, d_labels = create_pie_chart_dictionaries(catalogue_a, detailed=True)

    pie_chart(list(d.keys()), list(d.values()), f"{query_a.path}/out_stats/", f"{catalogue_a.id}_{query_a.id}_found_catalogue_genes_detailed", f"Detailed overlap scores of {catalogue_a.id} with {query_a.id} annotation", hovertext_labels=list(d_labels.values()), colours="extreme")

    if catalogue_expressed_path:
        location = f"{pickle_path}expressed_{catalogue_name}_annotation.pkl"
        if os.path.isfile(location):
            expressed_catalogue_a = pickle_load(location)
        else:
            expressed_catalogue_a = catalogue_pre.copy()
            expressed_catalogue_a.add_gene_names_and_descriptors(file_path=catalogue_expressed_path)
            expressed_catalogue_a.remove_noncoding_transcripts_from_coding_genes()
            expressed_catalogue_a.remove_genes_without_names()
            expressed_catalogue_a.id = f"expressed_{catalogue_name}"
            expressed_catalogue_a.name = expressed_catalogue_a.id.split("_on_")[0]
            expressed_catalogue_a.export_gff()
            expressed_catalogue_a.update_stats(export=True, genome=genomes[expressed_catalogue_a.genome])
            pickle_save(location, expressed_catalogue_a)

        expressed_catalogue_a.clear_overlaps()
        expressed_catalogue_a.detect_gene_overlaps(query_a)
        expressed_catalogue_a.export_equivalences(stringent=False, return_df=False, NAs=False, export_csv=True)

        d, d_labels = create_pie_chart_dictionaries(expressed_catalogue_a)

        pie_chart(list(d.keys()), list(d.values()), f"{query_a.path}/out_stats/", f"{expressed_catalogue_a.id}_{query_a.id}_found_catalogue_genes_simplified", f"Simplified overlap scores of {expressed_catalogue_a.id} with {query_a.id} annotation", hovertext_labels=list(d_labels.values()), colours="extreme")

        d, d_labels = create_pie_chart_dictionaries(expressed_catalogue_a, detailed=True)

        pie_chart(list(d.keys()), list(d.values()), f"{query_a.path}/out_stats/", f"{expressed_catalogue_a.id}_{query_a.id}_found_catalogue_genes_detailed", f"Detailed overlap scores of {expressed_catalogue_a.id} with {query_a.id} annotation", hovertext_labels=list(d_labels.values()), colours="extreme")

    if catalogue_pseudo_path:
        location = f"{pickle_path}pseudo_{catalogue_name}_annotation.pkl"
        if os.path.isfile(location):
            pseudo_catalogue_a = pickle_load(location)
        else:
            pseudo_catalogue_a = catalogue_pre.copy()
            pseudo_catalogue_a.add_gene_names_and_descriptors(file_path=catalogue_pseudo_path)
            pseudo_catalogue_a.remove_noncoding_transcripts_from_coding_genes()
            pseudo_catalogue_a.remove_genes_without_names()
            pseudo_catalogue_a.id = f"pseudo_{catalogue_name}"
            pseudo_catalogue_a.name = pseudo_catalogue_a.id.split("_on_")[0]
            pseudo_catalogue_a.export_gff()
            pseudo_catalogue_a.update_stats(export=True, genome=genomes[pseudo_catalogue_a.genome])
            pickle_save(location, pseudo_catalogue_a)

        pseudo_catalogue_a.clear_overlaps()
        pseudo_catalogue_a.detect_gene_overlaps(query_a)
        pseudo_catalogue_a.export_equivalences(stringent=False, return_df=False, NAs=False, export_csv=True)

        d, d_labels = create_pie_chart_dictionaries(pseudo_catalogue_a)

        pie_chart(list(d.keys()), list(d.values()), f"{query_a.path}/out_stats/", f"{pseudo_catalogue_a.id}_{query_a.id}_found_catalogue_genes_simplified", f"Simplified overlap scores of {pseudo_catalogue_a.id} with {query_a.id} annotation", hovertext_labels=list(d_labels.values()), colours="extreme")

        d, d_labels = create_pie_chart_dictionaries(pseudo_catalogue_a, detailed=True)

        pie_chart(list(d.keys()), list(d.values()), f"{query_a.path}/out_stats/", f"{pseudo_catalogue_a.id}_{query_a.id}_found_catalogue_genes_detailed", f"Detailed overlap scores of {pseudo_catalogue_a.id} with {query_a.id} annotation", hovertext_labels=list(d_labels.values()), colours="extreme")

# plot_catalogue_stats("final_reduced_on_40X_ref", "v3_from_12Xv2_raw_on_40X_ref", "catalogue_v3_on_40X_ref", "../../genomes_and_annotation/grapevine/PN40024/12Xv2_annotation/catalogue/vvi_lab_catalogue_V4_TM39_unofficial_alternative_stats.csv", "../../genomes_and_annotation/grapevine/PN40024/12Xv2_annotation/catalogue/vvi_lab_catalogue_V4_TM39_unofficial_alternative_expressed.csv", "../../genomes_and_annotation/grapevine/PN40024/12Xv2_annotation/catalogue/vvi_lab_catalogue_V4_TM39_unofficial_alternative_pseudo.csv")

# plot_catalogue_stats("EVM_PASA_on_40X_ref", "v3_from_12Xv2_raw_on_40X_ref", "catalogue_v3_on_40X_ref", "../../genomes_and_annotation/grapevine/PN40024/12Xv2_annotation/catalogue/vvi_lab_catalogue_V4_TM39_unofficial_alternative_stats.csv", "../../genomes_and_annotation/grapevine/PN40024/12Xv2_annotation/catalogue/vvi_lab_catalogue_V4_TM39_unofficial_alternative_expressed.csv", "../../genomes_and_annotation/grapevine/PN40024/12Xv2_annotation/catalogue/vvi_lab_catalogue_V4_TM39_unofficial_alternative_pseudo.csv")

# plot_catalogue_stats("final_reduced_without_PN40024_on_40X_ref", "v3_from_12Xv2_raw_on_40X_ref", "catalogue_v3_on_40X_ref", "../../genomes_and_annotation/grapevine/PN40024/12Xv2_annotation/catalogue/vvi_lab_catalogue_V4_TM39_unofficial_alternative_stats.csv", "../../genomes_and_annotation/grapevine/PN40024/12Xv2_annotation/catalogue/vvi_lab_catalogue_V4_TM39_unofficial_alternative_expressed.csv", "../../genomes_and_annotation/grapevine/PN40024/12Xv2_annotation/catalogue/vvi_lab_catalogue_V4_TM39_unofficial_alternative_pseudo.csv")

# plot_catalogue_stats("v4.1_on_40X_ref", "v3_from_12Xv2_raw_on_40X_ref", "catalogue_v3_on_40X_ref", "../../genomes_and_annotation/grapevine/PN40024/12Xv2_annotation/catalogue/vvi_lab_catalogue_V4_TM39_unofficial_alternative_stats.csv", "../../genomes_and_annotation/grapevine/PN40024/12Xv2_annotation/catalogue/vvi_lab_catalogue_V4_TM39_unofficial_alternative_expressed.csv", "../../genomes_and_annotation/grapevine/PN40024/12Xv2_annotation/catalogue/vvi_lab_catalogue_V4_TM39_unofficial_alternative_pseudo.csv")

# plot_catalogue_stats("v4.3_on_40X_ref", "v3_from_12Xv2_raw_on_40X_ref", "catalogue_v3_on_40X_ref", "../../genomes_and_annotation/grapevine/PN40024/12Xv2_annotation/catalogue/vvi_lab_catalogue_V4_TM39_unofficial_alternative_stats.csv", "../../genomes_and_annotation/grapevine/PN40024/12Xv2_annotation/catalogue/vvi_lab_catalogue_V4_TM39_unofficial_alternative_expressed.csv", "../../genomes_and_annotation/grapevine/PN40024/12Xv2_annotation/catalogue/vvi_lab_catalogue_V4_TM39_unofficial_alternative_pseudo.csv")

# plot_catalogue_stats("5.1_on_T2T_ref", "v3_from_12Xv2_raw_on_T2T_ref", "catalogue_v3_on_T2T_ref", "../../genomes_and_annotation/grapevine/PN40024/12Xv2_annotation/catalogue/vvi_lab_catalogue_V4_TM39_unofficial_alternative_stats.csv", "../../genomes_and_annotation/grapevine/PN40024/12Xv2_annotation/catalogue/vvi_lab_catalogue_V4_TM39_unofficial_alternative_expressed.csv", "../../genomes_and_annotation/grapevine/PN40024/12Xv2_annotation/catalogue/vvi_lab_catalogue_V4_TM39_unofficial_alternative_pseudo.csv")

plot_catalogue_stats("titan_on_malbec_prune_changed", "PN40024.v3_on_malbec_prune_changed", "catalogue_v3_on_malbec_prune_changed", "../../genomes_and_annotation/grapevine/PN40024/12Xv2_annotation/catalogue/vvi_lab_catalogue_V4_TM39_unofficial_alternative_stats.csv", "../../genomes_and_annotation/grapevine/PN40024/12Xv2_annotation/catalogue/vvi_lab_catalogue_V4_TM39_unofficial_alternative_expressed.csv", "../../genomes_and_annotation/grapevine/PN40024/12Xv2_annotation/catalogue/vvi_lab_catalogue_V4_TM39_unofficial_alternative_pseudo.csv")

now = time.time()
lapse = now - start
print(f"\nCreating and updating {species} annotation objects as well as exporting all features and gene lengths took {round((lapse/60)/60, 1)} hours\n")


