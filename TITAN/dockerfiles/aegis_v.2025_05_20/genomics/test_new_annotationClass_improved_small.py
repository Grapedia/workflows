#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import time
import os
from modules.tools import pickle_load, pickle_save
from modules.geneclasses import Annotation
from modules.geneclasses import Genome

# pn40024v4 = Genome('PN40024V4', '/media/tomslab2/Storage/Antonio/annotation_pipeline/data/assemblies/PN40024_40X_REF.chr_renamed_soft_masked.fasta')

# print('transcriptome evidences')
# stringtie_unstranded_default_STAR_2 = Annotation('stringtie_unstranded_default_STAR_2', '/media/tomslab2/Storage/Antonio/annotation_pipeline/Evidences_BLAST_hits/gffs/transcriptome_gffs/stringtie_unstranded_default_STAR_2.gff', pn40024v4)
# # stringtie_stranded_default_STAR_2 = Annotation('stringtie_stranded_default_STAR_2', '/media/tomslab2/Storage/Antonio/annotation_pipeline/Evidences_BLAST_hits/gffs/transcriptome_gffs/stringtie_stranded_default_STAR_2.gff', pn40024v4)

# print('merge transcripts evidence')
# trascript_evidence = stringtie_unstranded_default_STAR_2
# # trascript_evidence.merge_annotations(stringtie_stranded_default_STAR_2, overlap_threshold=0)

# trascript_evidence.generate_sequences(pn40024v4)
# trascript_evidence.rework_CDSs()
# trascript_evidence.generate_sequences(pn40024v4)

# print('export_gff merge')
# trascript_evidence.export_gff('/media/tomslab2/Storage/Antonio/annotation_pipeline/Evidences_BLAST_hits/gffs/test', 'trascript_evidence_test.gff3')
# trascript_evidence.export_proteins(only_main=True, verbose=False, custom_path='/media/tomslab2/Storage/Antonio/annotation_pipeline/Evidences_BLAST_hits/gffs/test')

#####
#
#   Test tomato overlaps 
#
#####

# tom_genome = Genome('tom_genome', '/media/tomslab2/Storage/Antonio/TOMSBioLab/genomes_and_annotation/tomato/Heinz_ITAG4_genome/ITAG4.fasta')
# tom_annot = Annotation('tom_annot', '/media/tomslab2/Storage/Antonio/TOMSBioLab/genomes_and_annotation/tomato/Heinz_ITAG4_annotation/ITAG4.0_priority_4.2_merge_tidy_fc_myb24_mod.gff3', tom_genome)
# tom_annot.detect_gene_overlaps()

# tom_annot.export_equivalences(export_self=True, custom_path='/media/tomslab2/Storage/Antonio/TOMSBioLab/repositorios/old_genomics/test_overlaps_new', export_csv=True, return_df=False, stringent=False, NAs=False)


#################
#
#
#    Test unstranded transcripts
#
#
################

pn40024v4 = Genome('PN40024V4', '/media/tomslab2/Storage/Antonio/annotation_pipeline/data/assemblies/PN40024_40X_REF.chr_renamed_soft_masked.fasta')

stringtie_stranded_default_STAR_2 = Annotation('stringtie_stranded_default_STAR_2', '/media/tomslab2/Storage/Antonio/annotation_pipeline/Evidences_BLAST_hits/gffs/transcriptome_gffs/stringtie_stranded_default_STAR_2.gff', pn40024v4)
stringtie_stranded_default_STAR_2.rework_CDSs(pn40024v4)

print()
