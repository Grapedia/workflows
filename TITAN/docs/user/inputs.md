# TITAN input preparation

## Table Of Contents

- [Required Inputs](#required-inputs)
- [Assemblies](#assemblies)
- [Previous Annotation](#previous-annotation)
- [RNA-seq Samplesheet](#rna-seq-samplesheet)
- [Protein Samplesheet](#protein-samplesheet)
- [EGAPx YAML](#egapx-yaml)
- [Input Validation](#input-validation)


This page documents the files required by TITAN and the expected formats. Use
absolute paths for production and HPC runs so paths remain valid on compute
nodes and inside container bind mounts.

## Required Inputs

| TITAN parameter | Required content |
| --- | --- |
| `--new_assembly` | Target genome assembly to annotate, FASTA. |
| `--previous_assembly` | Previous/reference assembly used by Liftoff, FASTA. |
| `--previous_annotations` | Annotation on the previous assembly, GFF3. |
| `--RNAseq_samplesheet` | CSV with RNA-seq inputs. |
| `--RNAseq_data_dir` | Directory containing local RNA-seq FASTQ/FASTA files referenced by the samplesheet. |
| `--protein_samplesheet` | CSV listing protein FASTA evidence files. |
| `--egapx_paramfile` | EGAPx YAML input file. |
| `--output_dir` | Output directory for TITAN results. |

Recommended layout:

```text
project/
  assemblies/
    previous.fa
    target.fa
  annotations/
    previous.gff3
  rnaseq/
    RNAseq_samplesheet.csv
    leaf_single.fastq.gz
    berry_paired_1.fastq.gz
    berry_paired_2.fastq.gz
    isoseq_leaf.fastq.gz
  proteins/
    protein_samplesheet.csv
    araport.fa
    swissprot_plants.fa
  egapx/
    input_egapx.yaml
  titan_out/
```

## Assemblies

`--new_assembly` is the target assembly TITAN annotates. `--previous_assembly`
is the assembly associated with `--previous_annotations` and is used by
Liftoff.

Both files should be valid FASTA files with stable sequence identifiers.
Sequence identifiers used by the previous GFF3 must match the previous assembly.
For production, keep assembly FASTA files immutable once a run has started.

## Previous Annotation

`--previous_annotations` must be a GFF3 annotation matching
`--previous_assembly`.

TITAN cleans Liftoff-derived GFF3 where needed before downstream AGAT/Aegis
steps, but input GFF3 should still be structurally valid and use sequence IDs
present in the previous assembly.

## RNA-seq Samplesheet

TITAN expects this CSV header:

```csv
sampleID,SRA_or_FASTQ,library_layout
```

Allowed `library_layout` values:

* `single`: single-end short reads.
* `paired`: paired-end short reads.
* `long`: long-read RNA-seq.

Allowed `SRA_or_FASTQ` values:

* `FASTQ`: local gzip-compressed FASTQ files under `RNAseq_data_dir`.
* `FASTA`: local FASTA file under `RNAseq_data_dir`; only valid with
  `library_layout=long`.
* `SRA`: SRR, ERR or DRR accession matching `sampleID`.

Example:

```csv
sampleID,SRA_or_FASTQ,library_layout
leaf_single,FASTQ,single
berry_paired,FASTQ,paired
isoseq_leaf,FASTQ,long
```

Local files are inferred from `sampleID` under `RNAseq_data_dir`:

```text
single: <sampleID>.fastq.gz
paired: <sampleID>_1.fastq.gz and <sampleID>_2.fastq.gz
long:   <sampleID>.fastq.gz or <sampleID>.fasta
```

SRA accessions are resolved through ENA FASTQ metadata. Production runs are
more reproducible when FASTQ/FASTA files are downloaded and staged outside
TITAN.

Long-read processing is automatic: if at least one row has
`library_layout=long`, TITAN launches the long-read branch and passes long-read
evidence to Aegis. FLAIR and SQANTI3 still require their own optional flags.

## Protein Samplesheet

TITAN expects this CSV header:

```csv
organism,filename
```

Example:

```csv
organism,filename
Araport,araport.fa
SwissProtPlants,swissprot_plants.fa
```

The `filename` entries should point to protein FASTA files available to
BRAKER3. Relative paths are resolved from the TITAN project directory.

Protein files must be multi-FASTA protein files suitable for ProtHint:

* every record has a header;
* every sequence is non-empty;
* the first word of each header is unique;
* terminal `*`, `.` and `-` are allowed and cleaned;
* internal `*` stop codons and non-protein characters fail validation.

## EGAPx YAML

EGAPx takes its own YAML input through `--egapx_paramfile`. TITAN validates
that the YAML contains `genome`, `taxid` and `organism`, and that `genome`
points to a readable FASTA file. Relative `genome` paths are resolved from the
YAML file directory.

Minimum example:

```yaml
genome: /absolute/path/to/project/assemblies/target.fa
taxid: 29760
organism: Vitis vinifera
```

Production runs should provide RNA-seq evidence to EGAPx as well:

```yaml
genome: /absolute/path/to/project/assemblies/target.fa
taxid: 29760
locus_tag_prefix: VITIS

short_reads:
  - - leaf_single
    - - /absolute/path/to/project/rnaseq/leaf_single.fastq.gz
  - - berry_paired
    - - /absolute/path/to/project/rnaseq/berry_paired_1.fastq.gz
      - /absolute/path/to/project/rnaseq/berry_paired_2.fastq.gz

long_reads:
  - - isoseq_leaf
    - - /absolute/path/to/project/rnaseq/isoseq_leaf.fastq.gz
```

Alternative for many short-read files:

```yaml
genome: /absolute/path/to/project/assemblies/target.fa
taxid: 29760
short_reads: /absolute/path/to/project/egapx/short_reads.txt
```

`short_reads.txt`:

```text
leaf_single /absolute/path/to/project/rnaseq/leaf_single.fastq.gz
berry_paired /absolute/path/to/project/rnaseq/berry_paired_1.fastq.gz
berry_paired /absolute/path/to/project/rnaseq/berry_paired_2.fastq.gz
```

Important EGAPx notes:

* `genome` should normally be an absolute path visible from the node running
  the EGAPx runner.
* `taxid` must be the NCBI taxonomy identifier for the target organism.
* Local RNA-seq files should be FASTA or FASTQ, not BAM.
* EGAPx performs its own masking; it does not require the EDTA-masked genome as
  input.
* TITAN launches EGAPx as a nested Nextflow workflow from the host environment,
  so the host must provide `python3`, `curl`, `tar` and the selected nested
  executor.

## Input Validation

The launcher runs `scripts/validate_inputs.py` before heavy execution. For
manual checks:

```bash
python3 scripts/validate_inputs.py \
  --project-dir "$PWD" \
  --new-assembly /absolute/path/to/project/assemblies/target.fa \
  --previous-assembly /absolute/path/to/project/assemblies/previous.fa \
  --previous-annotations /absolute/path/to/project/annotations/previous.gff3 \
  --rnaseq-samplesheet /absolute/path/to/project/rnaseq/RNAseq_samplesheet.csv \
  --rnaseq-data-dir /absolute/path/to/project/rnaseq \
  --protein-samplesheet /absolute/path/to/project/proteins/protein_samplesheet.csv \
  --egapx-paramfile /absolute/path/to/project/egapx/input_egapx.yaml \
  --egapx-executor singularity
```
