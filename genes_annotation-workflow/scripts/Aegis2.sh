#!/bin/bash

# Create diamond databases in case they are not created yet
module load diamond

diamond makedb --threads 20 --in /home/avelt/data2/2024_assembly_GW_RI_hifiasm/Riesling/2025_genes_annotation/workflows/genes_annotation-workflow/data/protein_data/arabidopsis_prot_2022_01.fasta \
-d /home/avelt/data2/2024_assembly_GW_RI_hifiasm/Riesling/2025_genes_annotation/workflows/genes_annotation-workflow/intermediate_files/BLAST_search_annotation_merge/arabidopsis_prot_2022_01.fasta

diamond makedb --threads 20 --in /home/avelt/data2/2024_assembly_GW_RI_hifiasm/Riesling/2025_genes_annotation/workflows/genes_annotation-workflow/data/protein_data/Viridiplantae_swissprot.fasta \
-d /home/avelt/data2/2024_assembly_GW_RI_hifiasm/Riesling/2025_genes_annotation/workflows/genes_annotation-workflow/intermediate_files/BLAST_search_annotation_merge/Viridiplantae_swissprot.fasta

diamond makedb --threads 20 --in /home/avelt/data2/2024_assembly_GW_RI_hifiasm/Riesling/2025_genes_annotation/workflows/genes_annotation-workflow/data/protein_data/eudicotyledons_uniprot.fasta \
-d /home/avelt/data2/2024_assembly_GW_RI_hifiasm/Riesling/2025_genes_annotation/workflows/genes_annotation-workflow/intermediate_files/BLAST_search_annotation_merge/eudicotyledons_uniprot.fasta


# Run diamond on the exported unique proteins

diamond blastp --threads 20 --db /home/avelt/data2/2024_assembly_GW_RI_hifiasm/Riesling/2025_genes_annotation/workflows/genes_annotation-workflow/intermediate_files/BLAST_search_annotation_merge/arabidopsis_prot_2022_01.fasta --ultra-sensitive --out "/home/avelt/data2/2024_assembly_GW_RI_hifiasm/Riesling/2025_genes_annotation/workflows/genes_annotation-workflow/intermediate_files/BLAST_search_annotation_merge/arabidopsis_prot_2022_01_vs_assembly.diamond" --outfmt 6 --query "/path/to/merge_annotation_on_PN40024_T2T_unique_proteins.fasta" --max-target-seqs 1 --evalue 1e-3
diamond blastp --threads 20 --db /home/avelt/data2/2024_assembly_GW_RI_hifiasm/Riesling/2025_genes_annotation/workflows/genes_annotation-workflow/intermediate_files/BLAST_search_annotation_merge/Viridiplantae_swissprot.fasta --ultra-sensitive --out "/home/avelt/data2/2024_assembly_GW_RI_hifiasm/Riesling/2025_genes_annotation/workflows/genes_annotation-workflow/intermediate_files/BLAST_search_annotation_merge/Viridiplantae_swissprot_vs_assembly.diamond" --outfmt 6 --query "/path/to/merge_annotation_on_PN40024_T2T_unique_proteins.fasta" --max-target-seqs 1 --evalue 1e-3
diamond blastp --threads 20 --db /home/avelt/data2/2024_assembly_GW_RI_hifiasm/Riesling/2025_genes_annotation/workflows/genes_annotation-workflow/intermediate_files/BLAST_search_annotation_merge/eudicotyledons_uniprot.fasta --ultra-sensitive --out "/home/avelt/data2/2024_assembly_GW_RI_hifiasm/Riesling/2025_genes_annotation/workflows/genes_annotation-workflow/intermediate_files/BLAST_search_annotation_merge/eudicotyledons_uniprot_vs_assembly.diamond" --outfmt 6 --query "/path/to/merge_annotation_on_PN40024_T2T_unique_proteins.fasta" --max-target-seqs 1 --evalue 1e-3

