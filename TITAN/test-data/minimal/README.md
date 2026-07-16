# TITAN minimal test data

Synthetic fixtures for fast local validation. These files are intentionally tiny and are not a biological benchmark for full genome annotation.

## Contents

Core inputs:

* `valid/reference.fa`: previous/reference assembly used by Liftoff.
* `valid/target.fa`: target assembly to annotate.
* `valid/reference.gff3`: previous annotation with gene, mRNA, exon and CDS features.
* `valid/input_egapx.yaml`: EGAPx parameter shape for config validation.
* `valid/rnaseq_samplesheet.csv`: RNA-seq samplesheet with `single`, `paired` and `long` layouts.
* `valid/rnaseq/*.fastq.gz`: tiny gzip-compressed FASTQ reads matching the RNA-seq samplesheet sample IDs.
* `valid/protein_samplesheet.csv`: protein evidence samplesheet with two organisms.
* `valid/proteins/*.fa`: minimal protein FASTA files used by BRAKER/Aegis fixtures.

Precomputed evidence fixtures:

* `valid/evidence/assembly_masked.EDTA.fasta`: EDTA-like masked assembly output.
* `valid/evidence/liftoff_previous_annotations.gff3`: Liftoff-like transferred annotation.
* `valid/evidence/augustus.hints.gff3`: BRAKER/AUGUSTUS-like GFF3.
* `valid/evidence/genemark.gtf`: GeneMark-like GTF.
* `valid/evidence/merged_star_stringtie_*.gtf`: STAR/StringTie short-read transcriptomes.
* `valid/evidence/merged_star_psiclass_*.gtf`: STAR/PsiCLASS transcriptomes.
* `valid/evidence/merged_minimap2_stringtie_long_reads_*.gtf`: Minimap2/StringTie long-read transcriptomes.

Negative fixtures:

* `invalid/empty.fa`: empty FASTA.
* `invalid/seqid_missing.gff3`: GFF3 seqid absent from FASTA.
* `invalid/bad_coordinates.gff3`: invalid coordinates.
* `invalid/bad_parent.gff3`: invalid Parent relation.

## Source and license

These fixtures are synthetic, created for TITAN development on 2026-07-16. They contain no production or confidential data.

## Expected use

Use these fixtures for static validation, parser tests, parameter checks and future `-stub-run` tests. They are too small for real EDTA, BRAKER3, STAR, HISAT2, Minimap2 or Aegis scientific validation.

Validate the dataset with:

```bash
python3 scripts/validate_minimal_test_data.py
sha256sum -c test-data/minimal/checksums.sha256
```
