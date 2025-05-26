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

# import _pydevd_bundle.pydevd_constants
# _pydevd_bundle.pydevd_constants.PYDEVD_WARN_EVALUATION_TIMEOUT = 30.

start = time.time()
os.system(f"mkdir -p {features_path}")
os.system(f"mkdir -p {pickle_path}")

genomes = {}

for k, v in genome_files.items():
    genomes[k] = Genome(k, v)
    genomes[k].export_feature_sizes()

annot_tags = ["genemark_on_40X_ref", "augustus_on_40X_ref", "v3_from_12Xv2_raw_on_40X_ref", "psiclass_stranded_on_40X_ref", "psiclass_unstranded_on_40X_ref", "stringtie_default_on_40X_ref", "stringtie_morus_on_40X_ref", "stringtie_stranded_default_on_40X_ref", "stringtie_stranded_morus_on_40X_ref", "stringtie_unstranded_on_40X_ref", "stringtie_unstranded_morus_on_40X_ref"]

annotations = {}

chosen_chromosome = "chr11"
chosen_coordinates = (7949000, 7963000)

for k, v in annotation_files.items():
    if k in annot_tags:
        location = f"{pickle_path}{chosen_chromosome}_{chosen_coordinates[0]}-{chosen_coordinates[1]}_{k}_annotation.pkl"
        if os.path.isfile(location):
            print(f"Loading {k} annotation object\n")
            annotations[k] = pickle_load(location)
        else:
            annot_name, genome = k.split("_on_")
            annotations[k] = Annotation(annot_name, v, genomes[genome], chosen_chromosome=chosen_chromosome, chosen_coordinates=chosen_coordinates)
            annotations[k].export_gff()
            pickle_save(location, annotations[k])

location = f"{pickle_path}{chosen_chromosome}_{chosen_coordinates[0]}-{chosen_coordinates[1]}_abinitio_transcriptome_merge_on_40X_ref_annotation.pkl"
if os.path.isfile(location):
    b = pickle_load(location)
else:
    b = annotations[annot_tags[0]].copy()
    del annotations[annot_tags[0]]

    for a in annotations.values():
        b.merge_annotations(a, overlap_threshold=0)

    b.name = f"{chosen_chromosome}_{chosen_coordinates[0]}-{chosen_coordinates[1]}_abinitio_transcriptome_merge"
    b.id = f"{chosen_chromosome}_{chosen_coordinates[0]}-{chosen_coordinates[1]}_abinitio_transcriptome_merge_on_40X_ref"
    b.calculate_transcript_masking(genomes[f"{b.genome}_hard"])
    b.export_gff(UTRs=True)
    pickle_save(location, b)

location = f"{pickle_path}full_{chosen_chromosome}_{chosen_coordinates[0]}-{chosen_coordinates[1]}_abinitio_transcriptome_merge_on_40X_ref_annotation.pkl"
if os.path.isfile(location):
    b = pickle_load(location)
else:
    b.name = f"full_{chosen_chromosome}_{chosen_coordinates[0]}-{chosen_coordinates[1]}_abinitio_transcriptome_merge"
    b.id = f"full_{chosen_chromosome}_{chosen_coordinates[0]}-{chosen_coordinates[1]}_abinitio_transcriptome_merge_on_40X_ref"
    b.make_alternative_transcripts_into_genes()
    b.export_unique_proteins(genome=genomes[b.genome])
    b.update_stats(export=True, genome=genomes[b.genome])
    b.generate_sequences(genome=genomes[b.genome])
    b.detect_gene_overlaps()
    # b.add_blast_hits("Araport", f"{b.path}/blast_results/David_test3_{chosen_chromosome}_{chosen_coordinates[0]}-{chosen_coordinates[1]}_Araport11.diamond")
    # b.add_blast_hits("Viridiplantae", f"{b.path}/blast_results/David_test3_{chosen_chromosome}_{chosen_coordinates[0]}-{chosen_coordinates[1]}_Viridiplantae.diamond")
    # b.add_blast_hits("Eudicotyledons", f"{b.path}/blast_results/David_test3_{chosen_chromosome}_{chosen_coordinates[0]}-{chosen_coordinates[1]}_Eudicotyledons.diamond")
    b.update()
    b.export_equivalences(stringent=False, return_df=False, export_csv=True, export_self=True, NAs=False)
    b.export_gff(UTRs=True)

    pickle_save(location, b)

location = f"{pickle_path}{chosen_chromosome}_{chosen_coordinates[0]}-{chosen_coordinates[1]}_abinitio_transcriptome_merge_reduced_on_40X_ref_annotation.pkl"
if os.path.isfile(location):
    b = pickle_load(location)
else:
    b.name = f"{chosen_chromosome}_{chosen_coordinates[0]}-{chosen_coordinates[1]}_abinitio_transcriptome_merge_reduced"
    b.id = f"{chosen_chromosome}_{chosen_coordinates[0]}-{chosen_coordinates[1]}_abinitio_transcriptome_merge_reduced_on_40X_ref"
    b.calculate_transcript_masking(genomes[f"{b.genome}_hard"])
    b.remove_redundancy(source_priority=["Araport", "Viridiplantae", "Eudicotyledons"], hard_masked_genome=genomes[f"{b.genome}_hard"])
    b.update_stats(export=True, genome=genomes[b.genome])
    b.export_gff(UTRs=False)
    b.clear_overlaps()
    b.detect_gene_overlaps()
    b.export_equivalences(stringent=False, return_df=False, export_csv=True, export_self=True, NAs=False)
    pickle_save(location, b)


location = f"{pickle_path}{chosen_chromosome}_{chosen_coordinates[0]}-{chosen_coordinates[1]}_abinitio_transcriptome_merge_reduced_combined_on_40X_ref_annotation.pkl"
if os.path.isfile(location):
    b = pickle_load(location)
else:
    b.name = f"{chosen_chromosome}_{chosen_coordinates[0]}-{chosen_coordinates[1]}_abinitio_transcriptome_merge_reduced_combined"
    b.id = f"{chosen_chromosome}_{chosen_coordinates[0]}-{chosen_coordinates[1]}_abinitio_transcriptome_merge_reduced_combined_on_40X_ref"
    b.combine_transcripts(genome=genomes[b.genome])
    b.update_stats(export=True, genome=genomes[b.genome])
    b.export_gff(UTRs=False)
    b.clear_overlaps()
    b.detect_gene_overlaps()
    b.export_equivalences(stringent=False, return_df=False, export_csv=True, export_self=True, NAs=False)
    pickle_save(location, b)

print('End')

# from modules.genefunctions import overlap

# print(b.chrs["chr01"]["STRG.951.1_955"].id, b.chrs["chr01"]["STRG.951.1_955"].start, b.chrs["chr01"]["STRG.951.1_955"].end, b.chrs["chr01"]["STRG.951.1_955"].size)
# print(b.chrs["chr01"]["chr01.6451.0"].id, b.chrs["chr01"]["chr01.6451.0"].start, b.chrs["chr01"]["chr01.6451.0"].end, b.chrs["chr01"]["chr01.6451.0"].size)
# print(b.chrs["chr01"]["chr01.6445.4"].id, b.chrs["chr01"]["chr01.6445.4"].start, b.chrs["chr01"]["chr01.6445.4"].end, b.chrs["chr01"]["chr01.6445.4"].size)

# b.chrs["chr01"]["STRG.951.1_955"].update()
# b.chrs["chr01"]["chr01.6451.0"].update()
# b.chrs["chr01"]["chr01.6445.4"].update()

# print(b.chrs["chr01"]["STRG.951.1_955"].id, b.chrs["chr01"]["STRG.951.1_955"].start, b.chrs["chr01"]["STRG.951.1_955"].end, b.chrs["chr01"]["STRG.951.1_955"].size)
# print(b.chrs["chr01"]["chr01.6451.0"].id, b.chrs["chr01"]["chr01.6451.0"].start, b.chrs["chr01"]["chr01.6451.0"].end, b.chrs["chr01"]["chr01.6451.0"].size)
# print(b.chrs["chr01"]["chr01.6445.4"].id, b.chrs["chr01"]["chr01.6445.4"].start, b.chrs["chr01"]["chr01.6445.4"].end, b.chrs["chr01"]["chr01.6445.4"].size)

# overlapping, overlap_bp = overlap(b.chrs["chr01"]["STRG.951.1_955"], b.chrs["chr01"]["chr01.6451.0"])
# print(overlapping, overlap_bp)
# overlapping, overlap_bp = overlap(b.chrs["chr01"]["STRG.951.1_955"], b.chrs["chr01"]["chr01.6445.4"]) 
# print(overlapping, overlap_bp)


# location = f"{pickle_path}chr01_abinitio_transcriptome_merge_reduced2_on_40X_ref_annotation.pkl"
# if os.path.isfile(location):
#     b = pickle_load(location)
# else:
#     b.remove_fully_intron_nested_genes()
#     b.select_best_possible_non_overlapping_UTR()
#     b.remove_transcript_redundancy()
#     b.remove_genes()
#     b.name = "chr01_abinitio_transcriptome_merge_reduced2"
#     b.id = "chr01_abinitio_transcriptome_merge_reduced2_on_40X_ref"
#     b.update()
#     b.export_gff(UTRs=True)
#     pickle_save(location, b)


# location = f"{pickle_path}chr01_abinitio_transcriptome_merge_reduced3_on_40X_ref_annotation.pkl"
# if os.path.isfile(location):
#     b = pickle_load(location)
# else:
#     b.correct_gene_transcript_and_subfeature_coordinates()
#     b.name = "chr01_abinitio_transcriptome_merge_reduced3"
#     b.id = "chr01_abinitio_transcriptome_merge_reduced3_on_40X_ref"
#     b.update()
#     b.export_gff(UTRs=True)
#     pickle_save(location, b)

now = time.time()
lapse = now - start
print(f"\nCreating and updating {species} annotation objects as well as other procedures took {round((lapse/60)/60, 1)} hours\n")