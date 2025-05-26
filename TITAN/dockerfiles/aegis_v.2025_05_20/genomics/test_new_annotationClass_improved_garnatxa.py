#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import time
import os
from modules.tools import pickle_load, pickle_save
from modules.geneclasses import Annotation
from modules.geneclasses import Genome

pn40024v4 = Genome('PN40024V4', '/storage/tom/Camille_T2T/transcriptomic_evidence/genome/PN40024_40X_REF.chr_renamed_soft_masked.fasta')

# print('ab initio evidences')
# augustus = Annotation('augustus', '/storage/tom/Camille_T2T/test_gffs/augustus.gff', pn40024v4)
# genemark = Annotation('genemark', '/storage/tom/Camille_T2T/test_gffs/genemark.gff', pn40024v4)
# # glimmerhmm = Annotation('glimmerhmm', '/storage/tom/Camille_T2T/test_gffs/glimmerhmm.gff', pn40024v4)
# # geneid = Annotation('geneid', '/storage/tom/Camille_T2T/test_gffs/geneid.gff', pn40024v4)
# liftoff = Annotation('liftoff', '/storage/tom/Camille_T2T/test_gffs/v3_copies_liftoff_ref_raw_filtered.gff3', pn40024v4)

# print('transcriptome evidences')
# psiclass_stranded_STAR = Annotation('psiclass_stranded_STAR', '/storage/tom/Camille_T2T/test_gffs/transcriptome_gffs/psiclass_stranded_STAR.gff', pn40024v4)
# psiclass_unstranded_STAR = Annotation('psiclass_unstranded_STAR', '/storage/tom/Camille_T2T/test_gffs/transcriptome_gffs/psiclass_unstranded_STAR.gff', pn40024v4)
# stringtie_stranded_default_STAR_2 = Annotation('stringtie_stranded_default_STAR_2', '/storage/tom/Camille_T2T/test_gffs/transcriptome_gffs/stringtie_stranded_default_STAR_2.gff', pn40024v4)
# stringtie_stranded_MorusCommands_STAR_2 = Annotation('stringtie_stranded_MorusCommands_STAR_2', '/storage/tom/Camille_T2T/test_gffs/transcriptome_gffs/stringtie_stranded_MorusCommands_STAR_2.gff', pn40024v4)
# stringtie_unstranded_default_STAR_2 = Annotation('stringtie_unstranded_default_STAR_2', '/storage/tom/Camille_T2T/test_gffs/transcriptome_gffs/stringtie_unstranded_default_STAR_2.gff', pn40024v4)
# stringtie_unstranded_MorusCommands_STAR_2 = Annotation('stringtie_unstranded_MorusCommands_STAR_2', '/storage/tom/Camille_T2T/test_gffs/transcriptome_gffs/stringtie_unstranded_MorusCommands_STAR_2.gff', pn40024v4)
# stringtie_SRA1_default_STAR_2 = Annotation('stringtie_SRA1_default_STAR_2', '/storage/tom/Camille_T2T/test_gffs/transcriptome_gffs/stringtie_SRA1_default_STAR_2.gff', pn40024v4)
# stringtie_SRA1_MorusCommands_STAR_2 = Annotation('stringtie_SRA1_MorusCommands_STAR_2', '/storage/tom/Camille_T2T/test_gffs/transcriptome_gffs/stringtie_SRA1_MorusCommands_STAR_2.gff', pn40024v4)

# print('merge transcripts evidence')
# trascript_evidence = psiclass_stranded_STAR.copy()
# trascript_evidence.merge_annotations(psiclass_unstranded_STAR, overlap_threshold=0)
# trascript_evidence.merge_annotations(stringtie_stranded_default_STAR_2, overlap_threshold=0)
# trascript_evidence.merge_annotations(stringtie_stranded_MorusCommands_STAR_2, overlap_threshold=0)
# trascript_evidence.merge_annotations(stringtie_unstranded_default_STAR_2, overlap_threshold=0)
# trascript_evidence.merge_annotations(stringtie_unstranded_MorusCommands_STAR_2, overlap_threshold=0)
# trascript_evidence.merge_annotations(stringtie_SRA1_default_STAR_2, overlap_threshold=0)
# trascript_evidence.merge_annotations(stringtie_SRA1_MorusCommands_STAR_2, overlap_threshold=0)

# print('merge all evidences')
# merge = augustus.copy()
# merge.merge_annotations(genemark, overlap_threshold=0)
# # merge.merge_annotations(glimmerhmm, overlap_threshold=0)
# # merge.merge_annotations(geneid, overlap_threshold=0)
# merge.merge_annotations(liftoff, overlap_threshold=0)
# merge.merge_annotations(trascript_evidence, overlap_threshold=0)

# merge.make_alternative_transcripts_into_genes()

# merge.id = 'merge_test11_on_40X_ref'
# merge.name = 'merge_test11'

# print('export_gff merge')
# merge.export_gff('/storage/tom/Camille_T2T/test_gffs/test', 'merge_test11.gff3')

# print('export_proteins merge')
# merge.export_unique_proteins(custom_path='/storage/tom/Camille_T2T/test_gffs/test', genome=pn40024v4)
# merge.update_stats(export=True, genome=pn40024v4)

# print('save pickle')
# pickle_save('/storage/tom/Camille_T2T/test_gffs/test/merge_test11.pkl', merge)




###################################################
#                                                 #
# Overlaps                                        #
#                                                 #
###################################################

merge = pickle_load('/storage/tom/Camille_T2T/test_gffs/test/merge_test11.pkl')

print('Masking...')
pn40024v4_hard_masked = Genome('pn40024v4_hard_masked', '/storage/tom/Camille_T2T/test_gffs/test/v4_genome_ref.fasta.mod.MAKER.masked')
merge.calculate_transcript_masking(hard_masked_genome = pn40024v4_hard_masked)

print('## Probably not needed in the future ## - Update objet to sort it')
merge.update()

print('Overlaps...')
merge.detect_gene_overlaps()

print('Adding BLAST results...')
merge.add_blast_hits('Araport', '/storage/tom/Camille_T2T/BLAST_search_annotation_merge/Araport11_pep_20220914_filt_merged_all11.diamond')
merge.add_blast_hits('Viridiplantae', '/storage/tom/Camille_T2T/BLAST_search_annotation_merge/Eudicotyledons_uniprot_2023.03.03_filt_merged_all11.diamond')
merge.add_blast_hits('Eudicots', '/storage/tom/Camille_T2T/BLAST_search_annotation_merge/Viridiplantae_swissprot_2023.03.02_filt_merged_all11.diamond')

print('Saving...')
pickle_save('/storage/tom/Camille_T2T/test_gffs/test/merge_test_overlaps11.pkl', merge)
merge.update()
merge.export_gff('/storage/tom/Camille_T2T/test_gffs/test', 'merge_test11.gff3')

print('Reduce redundancy...')
merge.remove_redundancy(source_priority = ['Araport', 'Viridiplantae', 'Eudicots'], hard_masked_genome=pn40024v4_hard_masked)

merge.id = 'final_annotation_test31_on_40X_ref'
merge.name = 'final_annotation_test31'

merge.export_gff('/storage/tom/Camille_T2T/test_gffs/test/', 'final_annotation_test31.gff3')
merge.export_equivalences(stringent=False, return_df=False, export_csv=True, export_self=True, NAs=False)
merge.generate_sequences(genome=pn40024v4, just_CDSs=True)
merge.export_proteins(only_main=True, verbose=False, custom_path="/storage/tom/Camille_T2T/test_gffs/test")
merge.export_proteins(only_main=False, verbose=False, custom_path="/storage/tom/Camille_T2T/test_gffs/test")