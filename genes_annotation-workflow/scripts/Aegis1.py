#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import time
import os
from modules.tools import pickle_load, pickle_save
from modules.annotation import Annotation
from modules.genome import Genome

# Load genome (first argument for a name of the genome and second argument the path)
PN40024_T2T = Genome('PN40024_T2T', '/storage/tom/PN40024_T2T.fasta')

# Define evidences names and paths (same as the genomes, a first argument for the name of the annotation and then the second for the path)
evidences = {
    'ab_initio': [
        ('augustus', '/storage/tom/augustus_T2T/Augustus/augustus.hints.gff3'),
        ('genemark', '/storage/tom/augustus_T2T/GeneMark-ETP/genemark.gff3'),
        ('liftoff', '/storage/tom/inigo/annotation/T2T/liftoff/soft_masked/T2T_ref_soft_masked.gff3')
    ],
    'transcriptome': [
        ('psiclass_stranded_T2T_STAR', '/storage/tom/T2T_mapping_bam/merged_T2T_PN40024_stranded_psiclass.gff3'),
        ('stringtie_stranded_default_T2T_STAR', '/storage/tom/T2T_mapping_bam/merged_T2T_PN40024_stranded_stringtie_default.gff3'),
        ('stringtie_stranded_AltCommands_T2T_STAR', '/storage/tom/T2T_mapping_bam/merged_T2T_PN40024_stranded_stringtie_morus.gff3'),
        ('psiclass_unstranded_T2T_STAR', '/storage/tom/T2T_mapping_bam/merged_T2T_PN40024_unstranded_psiclass.gff3'),
        ('stringtie_unstranded_default_T2T_STAR', '/storage/tom/T2T_mapping_bam/merged_T2T_PN40024_unstranded_stringtie_default.gff3'),
        ('stringtie_unstranded_AltCommands_T2T_STAR', '/storage/tom/T2T_mapping_bam/merged_T2T_PN40024_unstranded_stringtie_morus.gff3'),
        ('psiclass_SRA1_T2T_STAR', '/storage/tom/T2T_mapping_bam/SRA_T2T_transcriptome2_psiclass.gff3'),
        ('stringtie_SRA1_default_T2T_STAR', '/storage/tom/T2T_mapping_bam2/SRA_T2T_transcriptome2_stringtie_default.gff3'),
        ('stringtie_SRA1_AltCommands_T2T_STAR', '/storage/tom/T2T_mapping_bam2/SRA_T2T_transcriptome2_stringtie_morus.gff3'),
        ('stringtie_Isoseq_PN40024_default', '/storage/tom/T2T_mapping_bam/Isoseq/hq_transcripts.RI_rmv_on_T2T_sorted_stringtieDefaultCommands.gff3'),
        ('stringtie_Isoseq_PN40024_AltCommands', '/storage/tom/T2T_mapping_bam/Isoseq/hq_transcripts.RI_rmv_on_T2T_sorted_stringtieAltCommands.gff3')
    ]
}

# Create annotation objects
print('ab initio evidences')
ab_initio_annotations = [Annotation(name, path, PN40024_T2T) for name, path in evidences['ab_initio']]

# Create annotation objects
print('transcriptome evidences')
transcriptome_annotations = [Annotation(name, path, PN40024_T2T) for name, path in evidences['transcriptome']]

# Merge evidences
merged_annotation = transcriptome_annotations[0].copy()
for annotation in transcriptome_annotations[1:]:
    merged_annotation.merge(annotation)
    
for annotation in ab_initio_annotations:
    merged_annotation.merge(annotation)



merged_annotation = pickle_load('/storage/tom/v1_annotation_fixed.pkl')

merged_annotation.make_alternative_transcripts_into_genes()

merged_annotation.id = 'merge_annotation'
merged_annotation.name = 'merge_annotation_transcripts'

print('Exporting merged gff')
merged_annotation.export_gff('/storage/tom/', 'merged_annotation.gff3')

print('Exporting merged all unique proteins')
merged_annotation.export_unique_proteins(custom_path='/storage/tom/', genome=PN40024_T2T)
merged_annotation.update_stats(export=True, genome=PN40024_T2T)

# Save pickle (to load the merged object with all annotations in case anything fails in the following steps)
print('Save pickle')
pickle_save('/storage/tom/merged_annotation.pkl', merged_annotation)
