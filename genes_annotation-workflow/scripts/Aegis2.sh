#!/bin/bash

# Create diamond databases in case they are not created yet

# diamond makedb --threads 20 --in /storage/tom/BLAST_search_annotation_merge/Araport11_pep_20220914_filt.fasta -d /storage/tom/BLAST_search_annotation_merge
# diamond makedb --threads 20 --in /storage/tom/BLAST_search_annotation_merge/Viridiplantae_swissprot_2023.03.02_filt.fasta -d /storage/tom/BLAST_search_annotation_merge
# diamond makedb --threads 20 --in /storage/tom/BLAST_search_annotation_merge/Eudicotyledons_uniprot_2023.03.03_filt.fasta -d /storage/tom/BLAST_search_annotation_merge


# Run diamond on the exported unique proteins

diamond blastp --threads 20 --db /storage/tom/BLAST_search_annotation_merge/Araport11_pep_20220914_filt --ultra-sensitive --out "/storage/tom/BLAST_search_annotation_merge/Araport11_pep_20220914_filt_merged_all.diamond" --outfmt 6 --query "/storage/tom/merge_annotation_on_PN40024_T2T_unique_proteins.fasta" --max-target-seqs 1 --evalue 1e-3
diamond blastp --threads 20 --db /storage/tom/BLAST_search_annotation_merge/Viridiplantae_swissprot_2023.03.02_filt --ultra-sensitive --out "/storage/tom/BLAST_search_annotation_merge/Viridiplantae_swissprot_2023.03.02_filt_merged_all.diamond" --outfmt 6 --query "/storage/tom/merge_annotation_on_PN40024_T2T_unique_proteins.fasta" --max-target-seqs 1 --evalue 1e-3
diamond blastp --threads 20 --db /storage/tom/BLAST_search_annotation_merge/Eudicotyledons_uniprot_2023.03.03_filt --ultra-sensitive --out "/storage/tom/BLAST_search_annotation_merge/Eudicotyledons_uniprot_2023.03.03_filt_merged_all.diamond" --outfmt 6 --query "/storage/tom/merge_annotation_on_PN40024_T2T_unique_proteins.fasta" --max-target-seqs 1 --evalue 1e-3

