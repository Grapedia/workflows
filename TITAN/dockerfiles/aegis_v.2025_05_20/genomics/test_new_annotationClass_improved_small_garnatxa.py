#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import time
import os
from modules.tools import pickle_load, pickle_save
from modules.geneclasses import Annotation
from modules.geneclasses import Genome

pn40024v4 = Genome('PN40024V4', '/storage/tom/Camille_T2T/transcriptomic_evidence/genome/PN40024_40X_REF.chr_renamed_soft_masked.fasta')
pn40024v4_hard_masked = Genome('pn40024v4_hard_masked', '/storage/tom/Camille_T2T/transcriptomic_evidence/genome/PN40024_40X_REF.chr_renamed_hardmasked.fasta')
merge = pickle_load('/storage/tom/Camille_T2T/test_gffs/test/merge_test_overlaps12.pkl')

merge.temporary_add_missing_properties()
merge.update()

merge.calculate_transcript_masking(hard_masked_genome = pn40024v4_hard_masked)
merge.clear_overlaps()
merge.detect_gene_overlaps()


print('Reduce redundancy...')
merge.remove_redundancy(source_priority = ['Araport', 'Viridiplantae', 'Eudicots'], hard_masked_genome=pn40024v4_hard_masked)

merge.id = 'final_annotation_test42_on_40X_ref'
merge.name = 'final_annotation_test42'

merge.generate_sequences(genome=pn40024v4, just_CDSs=True)
merge.export_proteins(only_main=True, verbose=False, custom_path="/storage/tom/Camille_T2T/test_gffs/test")
merge.export_proteins(only_main=False, verbose=False, custom_path="/storage/tom/Camille_T2T/test_gffs/test")
merge.export_gff('/storage/tom/Camille_T2T/test_gffs/test/', 'final_annotation_test42.gff3')


# -----------------------------------

# pn40024v4 = Genome('PN40024V4', '/storage/tom/Camille_T2T/transcriptomic_evidence/genome/PN40024_40X_REF.chr_renamed_soft_masked.fasta')
# merge = Annotation('merge', '/storage/tom/Camille_T2T/test_gffs/test/final_annotation_test29.gff3', pn40024v4)
# v4_1 = Annotation('v4_1', '/storage/tom/Camille_T2T/test_gffs/test/v4_1_ref.gff3', pn40024v4)
# EvM = Annotation('EvM', '/storage/tom/Camille_T2T/test_gffs/test/EVM_unfiltered_V4_test.gff3', pn40024v4)

# merge.detect_gene_overlaps(other = v4_1)
# merge.export_equivalences(stringent=False, return_df=False, export_csv=True, export_self=False, NAs=False)
# os.rename("/storage/tom/Camille_T2T/test_gffs/test/overlaps/merge_on_PN40024V4_equivalences_noNAs.csv", "/storage/tom/Camille_T2T/test_gffs/test/overlaps/merge_overlap_v4.1.csv")
# merge.clear_overlaps()
# merge.detect_gene_overlaps(other = EvM)
# merge.export_equivalences(stringent=False, return_df=False, export_csv=True, export_self=False, NAs=False)
# os.rename("/storage/tom/Camille_T2T/test_gffs/test/overlaps/merge_on_PN40024V4_equivalences_noNAs.csv", "/storage/tom/Camille_T2T/test_gffs/test/overlaps/merge_overlap_EvM.csv")

# merge.id = "final_annotation_test29_on_40X_ref_combined"
# merge.name = "final_annotation_test29_combined"
# merge.combine_transcripts(genome=pn40024v4)
# merge.update()
# merge.export_gff()

# v4_1.id = "v4_1_on_40X_ref_combined"
# v4_1.name = "v4_1_combined"
# v4_1.combine_transcripts(genome=pn40024v4)
# v4_1.update()
# v4_1.export_gff()

# EvM.id = "EvM_on_40X_ref_combined"
# EvM.name = "EvM_combined"
# EvM.combine_transcripts(genome=pn40024v4)
# EvM.update()
# EvM.export_gff()