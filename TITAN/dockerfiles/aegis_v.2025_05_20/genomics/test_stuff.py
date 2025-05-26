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


# merge = pickle_load('/home/tomslab2/Descargas/test.pkl')
# for i, g in enumerate(merge.chrs['chr01'].values()):
#     if 'chr01.583.' in g.id:
#         print(i, g.id, g.start, g.end)

# merge.detect_gene_overlaps_test()
# print('end')

##########################

pn40024v4 = Genome('PN40024V4', '/media/tomslab2/Storage/Antonio/TOMSBioLab/genomes_and_annotation/grapevine/PN40024/40X_genome/v4_genome_ref.fasta')
# pn40024v4_hard_masked = Genome('pn40024v4_hard_masked', '/media/tomslab2/Storage/Antonio/annotation_pipeline/data/assemblies/PN40024_40X_REF.chr_renamed_hardmasked.fasta')

v4 = Annotation('v4-1', '/media/tomslab2/Storage/Antonio/TOMSBioLab/genomes_and_annotation/grapevine/PN40024/40X_annotation/v4_1_ref_tidy.gff3', pn40024v4)
# anno = Annotation('Amandine_pasa', '/media/tomslab2/Storage/Antonio/TOMSBioLab/genomes_and_annotation/grapevine/PN40024/40X_annotation_tests/Amandine_gff/Amandine_pasa.gff3', pn40024v4)
anno2 = Annotation('final_annotation_test41', '/home/tomslab2/Descargas/final_annotation_test41.gff3', pn40024v4)
# # anno.detect_gene_overlaps(v4)
anno2.detect_gene_overlaps(v4)
# # anno.export_equivalences(stringent=False, return_df=False, NAs=False, export_csv=True, custom_path="/media/tomslab2/Storage/Antonio/annotation_pipeline/eval/self")
anno2.export_equivalences(stringent=False, return_df=False, NAs=False, export_csv=True, custom_path="/media/tomslab2/Storage/Antonio/annotation_pipeline/eval/self")


# anno.combine_transcripts(pn40024v4)
# anno.export_gff(custom_path="/home/tomslab2/Descargas/", tag="Amandine_pasa_combined.gff3", main_only=False, UTRs=False)

# anno2.combine_transcripts(pn40024v4)
# anno2.export_gff(custom_path="/home/tomslab2/Descargas/", tag="final_annotation_test41_combined.gff3", main_only=False, UTRs=False)
