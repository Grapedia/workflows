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


pn40024v4 = Genome('PN40024V4', '/media/tomslab2/Storage/Antonio/TOMSBioLab/genomes_and_annotation/grapevine/PN40024/40X_genome/v4_genome_ref.fasta')
pn40024v4_hard_masked = Genome('pn40024v4_hard_masked', '/media/tomslab2/Storage/Antonio/annotation_pipeline/data/assemblies/PN40024_40X_REF.chr_renamed_hardmasked.fasta')

anno = Annotation('final_annotation_test36', '/media/tomslab2/Storage/Antonio/TOMSBioLab/genomes_and_annotation/grapevine/PN40024/40X_annotation_tests/Amandine_gff/Amandine_pasa.gff3', pn40024v4)

anno.calculate_transcript_masking(hard_masked_genome = pn40024v4_hard_masked)

CDS_masked_fraction_list = []
transcript_masked_fraction_list = []
for genes in anno.chrs.values():
    for g in genes.values():
        for t in g.transcripts.values():
            if t.main:
                transcript_masked_fraction = t.masked_fraction
                transcript_GC_content = t.gc_content
                for c in t.CDSs.values():
                    if c.main:
                        CDS_masked_fraction = c.masked_fraction
                        CDS_GC_content = c.gc_content
    
        CDS_masked_fraction_list.append(CDS_masked_fraction)
        transcript_masked_fraction_list.append(transcript_masked_fraction)

print(len([x for x in CDS_masked_fraction_list if x == 1]))
print(len(CDS_masked_fraction_list))
print((len([x for x in CDS_masked_fraction_list if x == 1]) / len(CDS_masked_fraction_list))*100, '%')


############## PieCharts ###############

def create_pie_chart_dictionaries(a):
    d = {}
    d["present"] = 0
    d["missing"] = 0

    d_labels = {}
    d_labels["present"] = []
    d_labels["missing"] = []

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


a = Annotation('v3', '/media/tomslab2/Storage/Antonio/TOMSBioLab/genomes_and_annotation/grapevine/PN40024/40X_annotation_transfer/original/v3_copies_liftoff_ref_raw.gff3', pn40024v4)
a.update_stats(export=False, genome=pn40024v4)

a1 = a.copy()
a1.add_gene_names_and_descriptors(file_path="/media/tomslab2/Storage/Antonio/TOMSBioLab/genomes_and_annotation/grapevine/PN40024/12Xv2_annotation/catalogue/vvi_lab_catalogue_V4_TM39_unofficial_alternative_stats.csv")
a1.remove_noncoding_transcripts_from_coding_genes()
a1.remove_genes_without_names()
a1.id = "catalogue_v3_on_40X_ref"
a1.name = "catalogue_v3"
a1.update_stats(export=False, genome=pn40024v4)


a_pseudo = a.copy()
a_pseudo.add_gene_names_and_descriptors(file_path="/media/tomslab2/Storage/Antonio/TOMSBioLab/genomes_and_annotation/grapevine/PN40024/12Xv2_annotation/catalogue/vvi_lab_catalogue_V4_TM39_unofficial_alternative_pseudo.csv")
a_pseudo.remove_noncoding_transcripts_from_coding_genes()
a_pseudo.remove_genes_without_names()
a_pseudo.id = "catalogue_v3_pseudo_on_40X_ref"
a_pseudo.name = "catalogue_v3_pseudo"
a_pseudo.update_stats(export=False, genome=pn40024v4)


a_expressed = a.copy()
a_expressed.add_gene_names_and_descriptors(file_path="/media/tomslab2/Storage/Antonio/TOMSBioLab/genomes_and_annotation/grapevine/PN40024/12Xv2_annotation/catalogue/vvi_lab_catalogue_V4_TM39_unofficial_alternative_expressed.csv")
a_expressed.remove_noncoding_transcripts_from_coding_genes()
a_expressed.remove_genes_without_names()
a_expressed.id = "catalogue_v3_expressed_on_40X_ref"
a_expressed.name = "catalogue_v3_expressed"
a_expressed.update_stats(export=False, genome=pn40024v4)

objects = [a, a1, a_expressed, a_pseudo]

anno.update_stats(export=True, genome=pn40024v4)

a1.detect_gene_overlaps(anno)
a1.export_equivalences(stringent=False, return_df=False, NAs=False, export_csv=True, custom_path="/media/tomslab2/Storage/Antonio/annotation_pipeline/eval/")
d, d_labels = create_pie_chart_dictionaries(a1)
pie_chart(list(d.keys()), list(d.values()), "/media/tomslab2/Storage/Antonio/annotation_pipeline/eval/", f"{a1.id}_{anno.id}_found_catalogue_genes_simplified", f"Simplified overlap scores of {a1.id} with {anno.id} annotation", hovertext_labels=list(d_labels.values()), colours="extreme")

a_pseudo.detect_gene_overlaps(anno)
a_pseudo.export_equivalences(stringent=False, return_df=False, NAs=False, export_csv=True, custom_path="/media/tomslab2/Storage/Antonio/annotation_pipeline/eval/")
d, d_labels = create_pie_chart_dictionaries(a_pseudo)
pie_chart(list(d.keys()), list(d.values()), "/media/tomslab2/Storage/Antonio/annotation_pipeline/eval/", f"{a_pseudo.id}_{anno.id}_found_catalogue_genes_simplified", f"Simplified overlap scores of {a_pseudo.id} with {anno.id} annotation", hovertext_labels=list(d_labels.values()), colours="extreme")

a_expressed.detect_gene_overlaps(anno)
a_expressed.export_equivalences(stringent=False, return_df=False, NAs=False, export_csv=True, custom_path="/media/tomslab2/Storage/Antonio/annotation_pipeline/eval/")
d, d_labels = create_pie_chart_dictionaries(a_expressed)
pie_chart(list(d.keys()), list(d.values()), "/media/tomslab2/Storage/Antonio/annotation_pipeline/eval/", f"{a_expressed.id}_{anno.id}_found_catalogue_genes_simplified", f"Simplified overlap scores of {a_expressed.id} with {anno.id} annotation", hovertext_labels=list(d_labels.values()), colours="extreme")

a1.clear_overlaps()
a_pseudo.clear_overlaps()
a_expressed.clear_overlaps()
