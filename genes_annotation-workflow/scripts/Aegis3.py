#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import time
import os
from modules.tools import pickle_load, pickle_save
from modules.annotation import Annotation
from modules.genome import Genome

merge = pickle_load('/storage/tom/merged_annotation.pkl')

print('Masking...')
pn40024T2T_hard_masked = Genome('pn40024T2T_hard_masked', '/storage/tom/T2T_genome_ref.fasta.mod.MAKER.masked')
merge.calculate_transcript_masking(hard_masked_genome = pn40024T2T_hard_masked)

print('## Probably not needed in the future ## - Update objet to sort it')
merge.update()

print('Overlaps...')
merge.detect_gene_overlaps()

print('Adding BLAST results...')
merge.add_blast_hits('Araport', '/storage/tom/BLAST_search_annotation_merge/Araport11_pep_20220914_filt_merged_all.diamond')
merge.add_blast_hits('Viridiplantae', '/storage/tom/BLAST_search_annotation_merge/Eudicotyledons_uniprot_2023.03.03_filt_merged_all.diamond')
merge.add_blast_hits('Eudicots', '/storage/tom/BLAST_search_annotation_merge/Viridiplantae_swissprot_2023.03.02_filt_merged_all.diamond')

# Save pickle (to load the object in case anything fails in the following steps)
print('Saving...')
pickle_save('/storage/tom/merged_annotation_blast.pkl', merge)

# Export gff to check BLAST results for each transcript
merge.update()
merge.export_gff('/storage/tom/inigo/araucaria/aegis/release', 'merged_annotation_blast.gff3')

print('Reduce redundancy...')
merge.remove_redundancy(source_priority = ['Araport', 'Viridiplantae', 'Eudicots'], hard_masked_genome=pn40024T2T_hard_masked)

# Save final pickle object
pickle_save('/storage/tom/final_annotation.pkl', merge)
merge.id = 'final_annotation'
merge.name = 'final_annotation'

merge.export_gff('/storage/tom/', 'final_annotation.gff3')
merge.export_equivalences(stringent=False, return_df=False, export_csv=True, export_self=True, NAs=False)
merge.generate_sequences(genome=PN40024_T2T, just_CDSs=True)
merge.export_proteins(only_main=True, verbose=False, custom_path="/storage/tom")
merge.export_proteins(only_main=False, verbose=False, custom_path="/storage/tom")
