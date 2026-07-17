# TITAN realistic lite fixture from PN40024 data

This fixture is generated from the real `data/` inputs, reduced to `chr01`
and a small but non-trivial RNA/protein evidence set.

Contents:

- target assembly: `assemblies/T2T_ref.chr01.fasta`
- previous assembly: `assemblies/v4_genome_ref.chr01.fasta`
- previous GFF3 annotation: `annotations/v4_3_just_ref.chr01.gff3`
- paired short reads: `depleted-stranded_Vitis14_S16_1/2.fastq.gz`
- single short reads: `PolyA_Vitis02_S1.fastq.gz`
- long reads: `PN40024_IsoSeq_chr01.fasta`
- reduced protein evidence: `5000` records from each cleaned source FASTA

Regenerate from the repository root with:

```bash
FASTQ_RECORDS=50000 LONG_RECORDS=3000 PROTEIN_RECORDS=5000 \
  scripts/create_realistic_lite_from_data.sh
```

Validate inputs with:

```bash
python3 scripts/validate_inputs.py \
  --project-dir . \
  --new-assembly data_example/realistic-lite/assemblies/T2T_ref.chr01.fasta \
  --previous-assembly data_example/realistic-lite/assemblies/v4_genome_ref.chr01.fasta \
  --previous-annotations data_example/realistic-lite/annotations/v4_3_just_ref.chr01.gff3 \
  --rnaseq-samplesheet data_example/realistic-lite/RNAseq_samplesheet.csv \
  --rnaseq-data-dir data_example/realistic-lite/rnaseq \
  --protein-samplesheet data_example/realistic-lite/protein_samplesheet.csv \
  --egapx-paramfile data_example/realistic-lite/egapx/input_egapx.yaml \
  --egapx-executor apptainer
```

Run a cluster smoke test first with `-stub-run`, then remove `-stub-run` for
the real end-to-end run:

```bash
nextflow -c data_example/realistic-lite/real_slurm_apptainer.config \
  run main.nf -profile slurm,apptainer -stub-run -ansi-log false
```
