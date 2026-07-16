# TITAN: The Intensive Transcript Annotation Pipeline

TITAN is a Nextflow pipeline for eukaryotic genome annotation. It combines transferred annotations, RNA-seq transcript evidence, protein-supported ab initio prediction, repeat masking, EGAPx annotation and final AEGIS integration.

Contributors: David Navarro, Antonio Santiago, Jose Tomas Matus, Amandine Velt, Camille Rustenholz and Marco Moretto.

Full setup details are in [docs/user/installation.md](docs/user/installation.md). Development audits and implementation notes are under [docs/development](docs/development).

## Current Contract

TITAN has one public execution mode: evidence generation, EDTA, EGAPx and AEGIS always run together in the same Nextflow graph. The former partial workflow modes are no longer supported, and `--workflow` is rejected at startup.

Mandatory steps:

* Liftoff transfers the previous annotation to the new assembly.
* EDTA creates the hard-masked target assembly required by AEGIS.
* EGAPx runs from `--egapx_paramfile` and contributes an additional GFF3 evidence.
* STAR/StringTie, STAR/PsiCLASS, HISAT2/StringTie and optional Minimap2/StringTie generate transcript evidence.
* BRAKER3 produces AUGUSTUS and GeneMark evidence.
* AEGIS merges named evidence channels and extracts final protein FASTA files.
* Diamond2GO annotates the final AEGIS protein sets.

Long-read processing is automatic when `RNAseq_samplesheet` contains at least one row with `library_layout=long`.

## Quick Start

From a fresh clone:

```bash
git clone git@github-amandine:Grapedia/workflows.git
cd workflows/TITAN
scripts/run-tests.sh
```

The quick test suite validates container pins, profile resolution, minimal fixtures, input validation and a full Nextflow `test` profile run in stub mode. It does not run scientific containers or validate biological annotation quality.

For production, prepare inputs and run through the launcher:

```bash
./launch_TITAN_example.sh \
  --profile slurm,apptainer \
  --output-dir /absolute/path/to/project/titan_out \
  --previous-assembly /absolute/path/to/project/assemblies/previous.fa \
  --new-assembly /absolute/path/to/project/assemblies/target.fa \
  --previous-annotations /absolute/path/to/project/annotations/previous.gff3 \
  --rnaseq-samplesheet /absolute/path/to/project/rnaseq/RNAseq_samplesheet.csv \
  --rnaseq-data-dir /absolute/path/to/project/rnaseq \
  --protein-samplesheet /absolute/path/to/project/proteins/protein_samplesheet.csv \
  --egapx-paramfile /absolute/path/to/project/egapx/input_egapx.yaml \
  --run-name titan_target_v1 \
  --resume
```

The launcher resolves paths, checks profile/container contracts, validates inputs and writes Nextflow reports under `${output_dir}/nextflow_reports`.

## Requirements

Use Linux or an HPC environment with:

* Bash, Git and Python 3.
* Java compatible with Nextflow.
* Nextflow `24.04.3` for the declared pipeline version.
* Docker for local container runs, or Apptainer/Singularity for HPC.
* `curl` and `tar` for the EGAPx runner bootstrap.

For EGAPx runner support, Python should provide `yaml`:

```bash
python3 -c 'import yaml; print("pyyaml OK")'
```

## Input Files

Required TITAN parameters:

| Parameter | Required content |
| --- | --- |
| `--new_assembly` | Target genome assembly to annotate, FASTA. |
| `--previous_assembly` | Reference assembly used by Liftoff, FASTA. |
| `--previous_annotations` | GFF3 annotation matching `--previous_assembly`. |
| `--RNAseq_samplesheet` | CSV with RNA-seq inputs. |
| `--RNAseq_data_dir` | Directory containing local FASTQ/FASTA files referenced by `sampleID`. |
| `--protein_samplesheet` | CSV listing protein FASTA evidence. |
| `--egapx_paramfile` | EGAPx YAML input file. |
| `--output_dir` | TITAN output directory. |

Recommended project layout:

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

TITAN validates FASTA, GFF3, samplesheets, local FASTQ/FASTA/protein file presence, EGAPx YAML basics, enum values and numeric options before heavy computation starts.

## RNA-Seq Samplesheet

Header:

```csv
sampleID,SRA_or_FASTQ,library_layout
```

Allowed `library_layout` values:

* `single`: single-end short reads.
* `paired`: paired-end short reads.
* `long`: long-read RNA-seq.

Allowed `SRA_or_FASTQ` values:

* `FASTQ`: local gzip-compressed FASTQ files under `RNAseq_data_dir`.
* `FASTA`: local FASTA file under `RNAseq_data_dir`, only valid with `library_layout=long`.
* `SRA`: accession handled by the historical preparation modules.

Example:

```csv
sampleID,SRA_or_FASTQ,library_layout
leaf_single,FASTQ,single
berry_paired,FASTQ,paired
isoseq_leaf,FASTQ,long
```

Local files are inferred from `sampleID`:

```text
single: <sampleID>.fastq.gz
paired: <sampleID>_1.fastq.gz and <sampleID>_2.fastq.gz
long:   <sampleID>.fastq.gz or <sampleID>.fasta
```

Production runs are more reproducible when FASTQ/FASTA files are downloaded and checked outside TITAN instead of relying on SRA expansion at runtime.

## Protein Samplesheet

Header:

```csv
organism,filename
```

Example:

```csv
organism,filename
Araport,/absolute/path/to/project/proteins/araport.fa
SwissProtPlants,/absolute/path/to/project/proteins/swissprot_plants.fa
```

Relative paths are resolved from the TITAN project directory by the validator. Absolute paths are recommended for production and HPC runs.

## EGAPx Input

`--egapx_paramfile` points to an EGAPx YAML file. TITAN validates at least `genome`, `taxid` and `organism`, and checks that `genome` points to a readable FASTA file.

Minimum:

```yaml
genome: /absolute/path/to/project/assemblies/target.fa
taxid: 29760
organism: Vitis vinifera
```

Recommended with RNA-seq evidence:

```yaml
genome: /absolute/path/to/project/assemblies/target.fa
taxid: 29760
organism: Vitis vinifera
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

EGAPx performs its own masking. It receives the original target assembly from its YAML, not the EDTA-masked assembly.

## Profiles

| Profile | Purpose |
| --- | --- |
| `test` | Local fixture-based validation; no Docker, Slurm or production data. |
| `local` | Local Docker-oriented production run. |
| `apptainer` | Apptainer/Singularity container runtime. |
| `slurm` | Slurm executor settings. |
| `slurm,apptainer` | Recommended HPC production profile combination. |
| `test,slurm,apptainer` | Configuration-resolution check using test data and HPC runtime settings. |

Useful checks:

```bash
nextflow config -profile test
nextflow config -profile local
nextflow config -profile slurm,apptainer
nextflow config -profile test,slurm,apptainer
```

Resource policy is centralized in [conf/base.config](conf/base.config). Active modules use process labels instead of local `cpus` directives:

| Label | Typical work |
| --- | --- |
| `process_low` | lightweight preparation, validation, Liftoff post-processing and provenance |
| `process_index` | genome or transcriptome index generation |
| `process_alignment` | read trimming and alignments |
| `process_transcriptome` | per-sample transcriptome assembly |
| `process_prediction` | EDTA, EGAPx and BRAKER3 prediction steps |
| `process_merge` | StringTie and GFFCompare merge steps |
| `process_aegis` | AEGIS and final functional annotation |

Override resources through `conf/base.config`, a profile config, or Nextflow process selectors. EDTA, EGAPx and Diamond2GO still honor `--edta_cpus`, `--egapx_cpus` and `--diamond2go_cpus` through `withName` selectors.

## Outputs

Main public output families:

| Output family | Main files | Location |
| --- | --- | --- |
| Liftoff | `liftoff_previous_annotations.gff3`, `unmapped_features.txt` | `${output_dir}` |
| EDTA | hard-masked assembly and TE files | `${output_dir}/tmp` |
| BRAKER3 | `augustus.hints.gff3`, `genemark.gtf`, `braker.gff3` | `${output_dir}` |
| Transcript evidence | merged STAR/StringTie, STAR/PsiCLASS and optional long-read GTFs | `${output_dir}` and intermediate transcriptome directories |
| EGAPx | `egapx.complete.genomic.gff3`, GTF, protein, CDS, transcript, ASN and `egapx_out/` | `${output_dir}/egapx` |
| AEGIS | `final_annotation.gff3`, `final_annotation_proteins_all.fasta`, `final_annotation_proteins_main.fasta` | `${output_dir}/aegis_outputs` |
| Diamond2GO | Diamond2GO annotation outputs | `${output_dir}/Diamond2GO_outputs` |
| Provenance | `evidence_manifest.json`, `versions.yml` | `${output_dir}/provenance` |
| Reports | Nextflow report, timeline, trace and DAG | `${output_dir}/nextflow_reports` when using the launcher |

`evidence_manifest.json` records main inputs, AEGIS evidence files, final AEGIS outputs, sizes and SHA-256 checksums. It is the current provenance record and the planned basis for more formal resume/reuse behavior in future work.

## Resume And Re-Runs

Use `-resume` or `./launch_TITAN_example.sh --resume` after interrupted runs. Keep the same `--output-dir`, `--work-dir`, profile combination and input paths when resuming.

The launcher protects non-empty output directories unless `--resume` or `--force` is provided. Use `--force` only when intentionally starting a new run into an existing directory.

For reproducible production runs, preserve:

* `${output_dir}/nextflow_reports`
* `${output_dir}/provenance/evidence_manifest.json`
* `${output_dir}/provenance/versions.yml`
* `.nextflow.log`
* the Nextflow work directory used for `-resume`

## Validation And CI

Local quick validation:

```bash
scripts/run-tests.sh
```

Targeted validation:

```bash
python3 scripts/validate_container_pins.py
python3 scripts/validate_profiles.py
python3 scripts/validate_minimal_test_data.py
python3 scripts/test_validate_inputs.py
nextflow run main.nf -profile test -stub-run -ansi-log false
```

GitHub Actions runs the same quick suite when files under `TITAN/**` or the TITAN CI workflow change.

## Troubleshooting

`Missing required parameter(s)`: provide every mandatory input and avoid empty shell variables in command lines.

`Required input file(s) not found`: verify paths from the machine or compute node running Nextflow. Prefer absolute paths for production.

`library_layout must be one of`: use only `single`, `paired` or `long` in `RNAseq_samplesheet`.

`at least one short-read row`: TITAN currently requires at least one `single` or `paired` short-read row.

`--workflow is no longer supported`: remove `--workflow`; TITAN always runs the full graph.

Nextflow version warning: TITAN declares `24.04.3`. Install or load that version for production if you want warning-free runs.

Apptainer or EGAPx image failures: confirm the pinned images in `nextflow.config` are reachable from compute nodes and that `TITAN_APPTAINER_CACHEDIR` points to writable shared storage.

STAR memory failures: increase `--STAR_memory_per_job` or adjust profile resources for high-depth RNA-seq.

Non-empty output directory: resume with `--resume`, or intentionally override with `--force`.

## Limitations

The `test` profile and CI run only stub execution on synthetic fixtures. They validate contracts, wiring and early failures, not scientific annotation quality.

Real Slurm and Apptainer behavior must still be validated on the target cluster because queue policy, shared filesystems and container cache settings are site-specific.

SRA inputs remain supported by historical modules, but local FASTQ/FASTA files are preferred for controlled production runs.

Automated final biological validation of the produced annotation is planned separately.

## Workflow Diagram

![Workflow Diagram](data_example/TITAN_diagram.jpg)

## Tool References

* [AEGIS](https://github.com/Tomsbiolab/aegis)
* [AGAT](https://github.com/NBISweden/AGAT)
* [BRAKER3](https://github.com/Gaius-Augustus/BRAKER)
* [Diamond2GO](https://github.com/rhysf/Diamond2GO)
* [EDTA](https://github.com/oushujun/EDTA)
* [EGAPx](https://github.com/ncbi/egapx)
* [fastp](https://github.com/OpenGene/fastp)
* [GFFCompare](https://ccb.jhu.edu/software/stringtie/gffcompare.shtml)
* [HISAT2](https://daehwankimlab.github.io/hisat2/)
* [Liftoff](https://github.com/agshumate/Liftoff)
* [Minimap2](https://github.com/lh3/minimap2)
* [PsiCLASS](https://github.com/splicebox/PsiCLASS)
* [Salmon](https://combine-lab.github.io/salmon/)
* [sra-tools](https://github.com/ncbi/sra-tools)
* [STAR](https://github.com/alexdobin/STAR)
* [StringTie](https://ccb.jhu.edu/software/stringtie/)
