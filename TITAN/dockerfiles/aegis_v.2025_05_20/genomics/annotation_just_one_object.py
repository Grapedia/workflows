#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import time
import os

from config.prunus import *
from modules.tools import pickle_load, pickle_save
from modules.annotation import Annotation
from modules.genome import Genome

aegis_output = f"{path}aegis_output/"
pickle_path = f"{aegis_output}pickles/"
features_path = f"{aegis_output}features/"
gff_path = f"{aegis_output}gffs/"
genome_path = f"{aegis_output}genomes/"
coordinates_path = f"{aegis_output}coordinates/"
stats_path = f"{aegis_output}stats/"

os.system(f"mkdir -p {features_path}")
os.system(f"mkdir -p {pickle_path}")
os.system(f"mkdir -p {gff_path}")
os.system(f"mkdir -p {genome_path}")
os.system(f"mkdir -p {coordinates_path}")
os.system(f"mkdir -p {stats_path}")

annotation_tag = "persica_on_persica"

genome_name = annotation_tag.split("_on_")[-1]
annotation_name = annotation_tag.split("_on_")[0]

g = Genome(genome_name, genome_files[genome_name], chromosome_dict=genome_dictionary[genome_name], rename_chromosomes=False)

location = f"{pickle_path}{annotation_tag}_annotation.pkl"
if os.path.isfile(location):
    print(f"Loading {annotation_tag} annotation object\n")
    a = pickle_load(location)
else:
    a = Annotation(annotation_name, annotation_files[annotation_tag], g)
    pickle_save(location, a)

a.generate_sequences(g)
a.export_proteins(unique_proteins_per_gene=False, custom_path=features_path, verbose=False, only_main=False)

print(f"{a.id} has {len(a.all_gene_ids.keys())} gene ids.")

print(f"{a.id} has {len(a.all_transcript_ids.keys())} transcript ids.")

main_transcript_count = 0

all_transcript_count = 0

coding_transcript_count = 0

CDS_count = 0

main_CDS_count = 0

CDS_seq_count = 0

for chrom, genes in a.chrs.items():
    for g in genes.values():
        for t in g.transcripts.values():
            if t.main:
                main_transcript_count += 1
            all_transcript_count += 1

            if t.coding:
                coding_transcript_count += 1

            for c in t.CDSs.values():
                if c.seq:
                    CDS_seq_count += 1
                if c.main:
                    main_CDS_count += 1
                CDS_count += 1


print(f"main_transcript_count={main_transcript_count}")

print(f"all_transcript_count={all_transcript_count}")

print(f"coding_transcript_count={coding_transcript_count}")

print(f"CDS_count={CDS_count}")

print(f"main_CDS_count={main_CDS_count}")

print(f"CDS_seq_count={CDS_seq_count}")