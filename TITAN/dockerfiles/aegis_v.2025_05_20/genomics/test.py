#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import time
import os
from modules.tools import pickle_load, pickle_save
from modules.geneclasses import Annotation
from modules.geneclasses import Genome, Annotation_geneless


pn40024v4 = Genome('PN40024V4', '/media/tomslab2/Storage/Antonio/annotation_pipeline/data/assemblies/PN40024_40X_REF.chr_renamed_soft_masked.fasta')
#abinitio_evidences = Annotation_geneless('abinitio_evidences', '/media/tomslab2/Storage/Antonio/annotation_pipeline/Evidences_BLAST_hits/ab_initio_evidences_corrected.gff', pn40024v4)

#pickle_save('/media/tomslab2/Storage/Antonio/annotation_pipeline/Evidences_BLAST_hits/ab_initio_evidences_corrected.pkl', abinitio_evidences)
#abinitio_evidences = pickle_load('/media/tomslab2/Storage/Antonio/annotation_pipeline/Evidences_BLAST_hits/ab_initio_evidences_corrected.pkl')
#abinitio_evidences.generate_sequences(pn40024v4)

#pickle_save('/media/tomslab2/Storage/Antonio/annotation_pipeline/Evidences_BLAST_hits/ab_initio_evidences_corrected_with_seqs.pkl', abinitio_evidences)

#abinitio_evidences = pickle_load('/media/tomslab2/Storage/Antonio/annotation_pipeline/Evidences_BLAST_hits/ab_initio_evidences_corrected_with_seqs.pkl')


#abinitio_evidences.export_proteins(False, True, custom_path='/media/tomslab2/Storage/Antonio/annotation_pipeline/Evidences_BLAST_hits/features_ab_initio_corrected')

#------------------------------

#transcript_evidence = Annotation_geneless('transcript_evidence', '/media/tomslab2/Storage/Antonio/annotation_pipeline/FINAL_OUTPUT_tmp/merged_transcriptomes.gff3', pn40024v4)
#transcript_evidence.generate_sequences(pn40024v4)
#pickle_save('/media/tomslab2/Storage/Antonio/annotation_pipeline/Evidences_BLAST_hits/merged_transcriptomes.pkl', transcript_evidence)
#transcript_evidence.generate_proteins_from_mRNA()
#transcript_evidence.export_transcript_proteins(False, True, custom_path='/media/tomslab2/Storage/Antonio/annotation_pipeline/Evidences_BLAST_hits/features_transcript_evidence')

# ----------------------------- Revisar el lunes
#abinitio_evidences = Annotation('abinitio_evidences', '/media/tomslab2/Storage/Antonio/annotation_pipeline/Evidences_BLAST_hits/gffs/ab_initio_evidences_corrected.gff', pn40024v4)
#abinitio_evidences.generate_sequences(pn40024v4)
#pickle_save('/media/tomslab2/Storage/Antonio/annotation_pipeline/Evidences_BLAST_hits/ab_initio_evidences_corrected.pkl', abinitio_evidences)

#abinitio_evidences.export_proteins(False, True, custom_path='/media/tomslab2/Storage/Antonio/annotation_pipeline/Evidences_BLAST_hits/features_ab_initio_corrected')

# -------------------------------
#pn40024v4 = Genome('PN40024V4', '/media/tomslab2/Storage/Antonio/annotation_pipeline/data/assemblies/v4_genome_merged.fasta')
#transcript_evidence = Annotation_geneless('transcript_evidence', '/media/tomslab2/Storage/Antonio/annotation_pipeline/test/test_merged_transcriptomes.gff3', pn40024v4)
#transcript_evidence.generate_sequences(pn40024v4)
#transcript_evidence.generate_proteins_from_mRNA()
#transcript_evidence.export_gff('/media/tomslab2/Storage/Antonio/annotation_pipeline/Evidences_BLAST_hits/gffs', 'transcript_evidences_CDSannotated.gff3')

# ------------------------------

#abinitio_evidences = pickle_load('/media/tomslab2/Storage/Antonio/annotation_pipeline/Evidences_BLAST_hits/ab_initio_evidences_corrected.pkl')
#transcript_evidences = pickle_load('/media/tomslab2/Storage/Antonio/annotation_pipeline/Evidences_BLAST_hits/gffs/transcript_evidences_CDSannotated_corrected.pkl')

#print('1')
#abinitio_evidences.merge_annotations(transcript_evidences, overlap_threshold=0)
#print('2')
#abinitio_evidences.generate_sequences(pn40024v4)
#print('3')
#pickle_save('/media/tomslab2/Storage/Antonio/annotation_pipeline/Evidences_BLAST_hits/gffs/merge_all.pkl', abinitio_evidences)
#print('4')
#abinitio_evidences.export_proteins(only_main=True, verbose=False, custom_path='/media/tomslab2/Storage/Antonio/annotation_pipeline/Evidences_BLAST_hits')
#print('5')
#abinitio_evidences.export_gff(custom_path='/media/tomslab2/Storage/Antonio/annotation_pipeline/Evidences_BLAST_hits/gffs')

# ----------------------------

#abinitio_evidences = pickle_load('/media/tomslab2/Storage/Antonio/annotation_pipeline/Evidences_BLAST_hits/gffs/merge_all.pkl')
#abinitio_evidences.detect_gene_overlaps_huge_file()
#pickle_save('/media/tomslab2/Storage/Antonio/annotation_pipeline/Evidences_BLAST_hits/gffs/merge_all_overlaps.pkl', abinitio_evidences)

# -----------------------------

abinitio_evidences = pickle_load('/media/tomslab2/Storage/Antonio/annotation_pipeline/Evidences_BLAST_hits/gffs/merge_all_overlaps.pkl')

print('Add Blast hits')
abinitio_evidences.add_blast_hits('Araport', '/media/tomslab2/Storage/Antonio/annotation_pipeline/Evidences_BLAST_hits/output_diamond/Araport11_pep_20220914_filt_merged_all.diamond', '/media/tomslab2/Storage/Antonio/annotation_pipeline/Evidences_BLAST_hits/repeated_sequences/repeated_sequences_merged_all.tsv')
abinitio_evidences.add_blast_hits('Viridiplantae', '/media/tomslab2/Storage/Antonio/annotation_pipeline/Evidences_BLAST_hits/output_diamond/Viridiplantae_swissprot_2023.03.02_filt_merged_all.diamond', '/media/tomslab2/Storage/Antonio/annotation_pipeline/Evidences_BLAST_hits/repeated_sequences/repeated_sequences_merged_all.tsv')
abinitio_evidences.add_blast_hits('Eudicots', '/media/tomslab2/Storage/Antonio/annotation_pipeline/Evidences_BLAST_hits/output_diamond/Eudicotyledons_uniprot_2023.03.03_filt_merged_all.diamond', '/media/tomslab2/Storage/Antonio/annotation_pipeline/Evidences_BLAST_hits/repeated_sequences/repeated_sequences_merged_all.tsv')

print('Reduce Redundancy')
abinitio_evidences.remove_redundancy(source_priority=['Araport', 'Viridiplantae', 'Eudicots'])
abinitio_evidences.id = 'final_annotation_reduced_redundancy_on_40X_test10'
abinitio_evidences.name = 'final_annotation_reduced_redundancy_test10'
pickle_save('/media/tomslab2/Storage/Antonio/annotation_pipeline/Evidences_BLAST_hits/gffs/final_annotation_reduced_redundancy.pkl', abinitio_evidences)

print('Export')
abinitio_evidences.export_gff(custom_path='/media/tomslab2/Storage/Antonio/annotation_pipeline/Evidences_BLAST_hits/gffs/preremove')
abinitio_evidences.remove_missing_transcripts_parent_references()
abinitio_evidences.export_gff(custom_path='/media/tomslab2/Storage/Antonio/annotation_pipeline/Evidences_BLAST_hits/gffs/postremove')

# ----------------------------

# abinitio_evidences = pickle_load('/media/tomslab2/Storage/Antonio/annotation_pipeline/Evidences_BLAST_hits/gffs/merge_all_overlaps.pkl')
# abinitio_evidences.id = 'pre_final_annotation_on_40X'
# abinitio_evidences.name = 'pre_final_annotation'
# abinitio_evidences.export_proteins(only_main=True, verbose=False, custom_path='/media/tomslab2/Storage/Antonio/annotation_pipeline/Evidences_BLAST_hits')
# abinitio_evidences.remove_missing_transcripts_parent_references()
# abinitio_evidences.export_gff(custom_path='/media/tomslab2/Storage/Antonio/annotation_pipeline/Evidences_BLAST_hits/gffs')


# abinitio_evidences = pickle_load('/media/tomslab2/Storage/Antonio/annotation_pipeline/Evidences_BLAST_hits/gffs/final_annotation_reduced_redundancy.pkl')
abinitio_evidences.generate_sequences(pn40024v4)
abinitio_evidences.export_proteins(only_main=True, verbose=False, custom_path='/media/tomslab2/Storage/Antonio/annotation_pipeline/Evidences_BLAST_hits')
