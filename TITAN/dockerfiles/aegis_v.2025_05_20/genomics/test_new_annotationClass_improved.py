#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import time
import os
from modules.tools import pickle_load, pickle_save
from modules.geneclasses import Annotation
from modules.geneclasses import Genome

# pn40024v4 = Genome('PN40024V4', '/media/tomslab2/Storage/Antonio/annotation_pipeline/data/assemblies/PN40024_40X_REF.chr_renamed_soft_masked.fasta')

# print('ab initio evidences')
# augustus = Annotation('augustus', '/media/tomslab2/Storage/Antonio/annotation_pipeline/Evidences_BLAST_hits/gffs/augustus.gff', pn40024v4)
# glimmerhmm = Annotation('glimmerhmm', '/media/tomslab2/Storage/Antonio/annotation_pipeline/Evidences_BLAST_hits/gffs/glimmerhmm.gff', pn40024v4)
# geneid = Annotation('geneid', '/media/tomslab2/Storage/Antonio/annotation_pipeline/Evidences_BLAST_hits/gffs/geneid.gff', pn40024v4)
# liftoff = Annotation('liftoff', '/media/tomslab2/Storage/Antonio/annotation_pipeline/Evidences_BLAST_hits/gffs/v3_copies_liftoff_ref_raw_filtered.gff3', pn40024v4)

# print('transcriptome evidences')
# psiclass_stranded_STAR = Annotation('psiclass_stranded_STAR', '/media/tomslab2/Storage/Antonio/annotation_pipeline/Evidences_BLAST_hits/gffs/transcriptome_gffs/psiclass_stranded_STAR.gff', pn40024v4)
# psiclass_unstranded_STAR = Annotation('psiclass_unstranded_STAR', '/media/tomslab2/Storage/Antonio/annotation_pipeline/Evidences_BLAST_hits/gffs/transcriptome_gffs/psiclass_unstranded_STAR.gff', pn40024v4)
# stringtie_stranded_default_STAR_2 = Annotation('stringtie_stranded_default_STAR_2', '/media/tomslab2/Storage/Antonio/annotation_pipeline/Evidences_BLAST_hits/gffs/transcriptome_gffs/stringtie_stranded_default_STAR_2.gff', pn40024v4)
# stringtie_stranded_MorusCommands_STAR_2 = Annotation('stringtie_stranded_MorusCommands_STAR_2', '/media/tomslab2/Storage/Antonio/annotation_pipeline/Evidences_BLAST_hits/gffs/transcriptome_gffs/stringtie_stranded_MorusCommands_STAR_2.gff', pn40024v4)
# stringtie_unstranded_default_STAR_2 = Annotation('stringtie_unstranded_default_STAR_2', '/media/tomslab2/Storage/Antonio/annotation_pipeline/Evidences_BLAST_hits/gffs/transcriptome_gffs/stringtie_unstranded_default_STAR_2.gff', pn40024v4)
# stringtie_unstranded_MorusCommands_STAR_2 = Annotation('stringtie_unstranded_MorusCommands_STAR_2', '/media/tomslab2/Storage/Antonio/annotation_pipeline/Evidences_BLAST_hits/gffs/transcriptome_gffs/stringtie_unstranded_MorusCommands_STAR_2.gff', pn40024v4)
# stringtie_SRA1_default_STAR_2 = Annotation('stringtie_SRA1_default_STAR_2', '/media/tomslab2/Storage/Antonio/annotation_pipeline/Evidences_BLAST_hits/gffs/transcriptome_gffs/stringtie_SRA1_default_STAR_2.gff', pn40024v4)
# stringtie_SRA1_MorusCommands_STAR_2 = Annotation('stringtie_SRA1_MorusCommands_STAR_2', '/media/tomslab2/Storage/Antonio/annotation_pipeline/Evidences_BLAST_hits/gffs/transcriptome_gffs/stringtie_SRA1_MorusCommands_STAR_2.gff', pn40024v4)


# print('merge transcripts evidence')
# trascript_evidence = psiclass_stranded_STAR.copy()
# del psiclass_stranded_STAR
# trascript_evidence.merge_annotations(psiclass_unstranded_STAR, overlap_threshold=0)
# del psiclass_unstranded_STAR
# trascript_evidence.merge_annotations(stringtie_stranded_default_STAR_2, overlap_threshold=0)
# del stringtie_stranded_default_STAR_2
# trascript_evidence.merge_annotations(stringtie_stranded_MorusCommands_STAR_2, overlap_threshold=0)
# del stringtie_stranded_MorusCommands_STAR_2
# trascript_evidence.merge_annotations(stringtie_unstranded_default_STAR_2, overlap_threshold=0)
# del stringtie_unstranded_default_STAR_2
# trascript_evidence.merge_annotations(stringtie_unstranded_MorusCommands_STAR_2, overlap_threshold=0)
# del stringtie_unstranded_MorusCommands_STAR_2
# trascript_evidence.merge_annotations(stringtie_SRA1_default_STAR_2, overlap_threshold=0)
# del stringtie_SRA1_default_STAR_2
# trascript_evidence.merge_annotations(stringtie_SRA1_MorusCommands_STAR_2, overlap_threshold=0)
# del stringtie_SRA1_MorusCommands_STAR_2
# print('rework CDSs on merge')
# trascript_evidence.rework_CDSs(pn40024v4)


# print('merge all evidences')
# merge = augustus.copy()
# del augustus
# merge.merge_annotations(glimmerhmm, overlap_threshold=0)
# del glimmerhmm
# merge.merge_annotations(geneid, overlap_threshold=0)
# del geneid
# merge.merge_annotations(liftoff, overlap_threshold=0)
# del liftoff
# merge.merge_annotations(trascript_evidence, overlap_threshold=0)
# del trascript_evidence

# merge.id = 'merge_test_on_40X_ref'
# merge.name = 'merge_test'

# merge.generate_sequences(pn40024v4, just_CDSs=True)

# print('export_gff merge')
# merge.export_gff('/media/tomslab2/Storage/Antonio/annotation_pipeline/Evidences_BLAST_hits/gffs/test', 'merge_test.gff3')

# print('export_proteins merge')
# merge.export_unique_proteins(custom_path='/media/tomslab2/Storage/Antonio/annotation_pipeline/Evidences_BLAST_hits/gffs/test')

# print('save pickle')
# pickle_save('/media/tomslab2/Storage/Antonio/annotation_pipeline/Evidences_BLAST_hits/new_pkls/merge_test.pkl', merge)


###################################################
#                                                 #
# Overlaps                                        #
#                                                 #
###################################################

# merge = pickle_load('/media/tomslab2/Storage/Antonio/annotation_pipeline/Evidences_BLAST_hits/new_pkls/merge_test.pkl')

# merge.clear_sequences()

merge = pickle_load('/media/tomslab2/Storage/Antonio/annotation_pipeline/Evidences_BLAST_hits/gffs/test/merge_test_no_sequences.pkl')

merge.detect_gene_overlaps()

merge.add_blast_hits('Araport', '/media/tomslab2/Storage/Antonio/annotation_pipeline/Evidences_BLAST_hits/BLAST/Araport11_pep_20220914_filt_merged_all.diamond')
merge.add_blast_hits('Viridiplantae', '/media/tomslab2/Storage/Antonio/annotation_pipeline/Evidences_BLAST_hits/BLAST/Eudicotyledons_uniprot_2023.03.03_filt.diamond')
merge.add_blast_hits('Eudicots', '/media/tomslab2/Storage/Antonio/annotation_pipeline/Evidences_BLAST_hits/BLAST/Viridiplantae_swissprot_2023.03.02_filt_merged_all.diamond')

pickle_save('/media/tomslab2/Storage/Antonio/annotation_pipeline/Evidences_BLAST_hits/new_pkls/bullshit.pkl', merge)

merge.remove_redundancy(source_priority = ['Araport', 'Viridiplantae', 'Eudicots'])
merge.export_gff('/media/tomslab2/Storage/Antonio/annotation_pipeline/Evidences_BLAST_hits/gffs/test', 'final_annotation_test11.gff3')
# pickle_save('/media/tomslab2/Storage/Antonio/annotation_pipeline/Evidences_BLAST_hits/new_pkls/merge_test_overlaps.pkl', merge)