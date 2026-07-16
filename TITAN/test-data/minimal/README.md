# TITAN minimal test data

Synthetic fixtures for fast local validation. These files are intentionally tiny and are not a biological benchmark for full genome annotation.

## Contents

* `valid/reference.fa`: previous/reference assembly.
* `valid/target.fa`: target assembly to annotate.
* `valid/reference.gff3`: minimal valid GFF3 with gene, mRNA, exon and CDS features.
* `valid/proteins.fa`: minimal protein FASTA.
* `valid/rnaseq_samplesheet.csv`: minimal samplesheet shape for TITAN parsing.
* `valid/protein_samplesheet.csv`: minimal protein samplesheet.
* `invalid/empty.fa`: empty FASTA.
* `invalid/seqid_missing.gff3`: GFF3 seqid absent from FASTA.
* `invalid/bad_coordinates.gff3`: invalid coordinates.
* `invalid/bad_parent.gff3`: invalid Parent relation.

## Source and license

These fixtures are synthetic, created for TITAN development on 2026-07-16. They contain no production or confidential data.

## Expected use

Use these fixtures for static validation, parser tests and future `-stub-run` tests. They are too small for real EDTA, BRAKER3, STAR, HISAT2, Minimap2 or Aegis scientific validation.
