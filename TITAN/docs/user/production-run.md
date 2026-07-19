# TITAN production run

This page documents the recommended production launch path, Slurm/Apptainer
usage, resume behavior and final checklist. Install TITAN first with
[installation.md](installation.md), then prepare inputs and reference data with
[inputs.md](inputs.md) and [reference-data.md](reference-data.md).

## Launcher

Use `launch_TITAN_example.sh` for production-oriented runs. It validates inputs,
checks profile/container contracts, resolves paths to absolute paths, prepares
supported reference data when requested and then runs Nextflow.

Base command:

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

The base command runs the mandatory graph: Liftoff, RNA-seq evidence, EDTA,
EGAPx, BRAKER3, AEGIS, Diamond2GO and default quality checks.

## Optional Branches

The launcher can both prepare and enable these branches:

```bash
./launch_TITAN_example.sh \
  ... \
  --prepare-egapx-cache \
  --enable-eggnog-mapper --prepare-eggnog-data \
  --enable-helixer --prepare-helixer-model --helixer-lineage land_plant \
  --enable-interproscan --prepare-interproscan-data \
  --enable-rfam --prepare-rfam-data \
  --enable-omark --prepare-omark-data \
  --enable-lncrna --prepare-cpat-model \
  --resume
```

Options without a dedicated launcher flag can be passed after `--` as native
Nextflow parameters:

```bash
./launch_TITAN_example.sh ... -- \
  --run_trnascan true \
  --run_busco true \
  --busco_data_dir /absolute/path/to/project/busco_data \
  --busco_lineage eudicotyledons_odb12.2 \
  --run_mikado true \
  --run_flair true \
  --run_sqanti3 true
```

See [reference-data.md](reference-data.md) before enabling branches that need
external databases.

## Profiles

Common profile checks:

```bash
nextflow config -profile local
nextflow config -profile apptainer
nextflow config -profile slurm,apptainer
nextflow config -profile test,slurm,apptainer
```

Use `test,slurm,apptainer` only for configuration resolution on minimal
fixtures. Production HPC runs should use `slurm,apptainer`.

Runtime resources are centralized in `conf/base.config` through process labels
such as `process_index`, `process_alignment`, `process_transcriptome`,
`process_prediction`, `process_merge` and `process_aegis`. Adjust resources in
config profiles or with Nextflow process selectors rather than editing module
files.

## Slurm and Apptainer

For Apptainer runs, keep `TITAN_APPTAINER_CACHEDIR` on a writable shared
filesystem visible to all compute nodes. The launcher sets it automatically to
`<output-dir>/apptainer-cache` when the selected profile includes `apptainer`,
unless it is already defined.

EGAPx is a nested Nextflow workflow. The `apptainer` profile sets
`egapx_executor = singularity` for the nested EGAPx run. On HPC without Docker,
also pass a prepared `--egapx_config_dir` if your site requires specific nested
EGAPx executor settings.

Site-specific examples are kept under [examples/colmar](../../examples/colmar).
Treat those as templates, not generic defaults.

## Direct Nextflow Command

The launcher is recommended, but the equivalent direct command is:

```bash
nextflow run main.nf \
  -profile slurm,apptainer \
  -name titan_target_v1 \
  -work-dir /absolute/path/to/project/titan_out/work \
  -with-dag /absolute/path/to/project/titan_out/nextflow_reports/titan_target_v1.dag.html \
  -ansi-log false \
  --output_dir /absolute/path/to/project/titan_out \
  --previous_assembly /absolute/path/to/project/assemblies/previous.fa \
  --new_assembly /absolute/path/to/project/assemblies/target.fa \
  --previous_annotations /absolute/path/to/project/annotations/previous.gff3 \
  --RNAseq_samplesheet /absolute/path/to/project/rnaseq/RNAseq_samplesheet.csv \
  --RNAseq_data_dir /absolute/path/to/project/rnaseq \
  --protein_samplesheet /absolute/path/to/project/proteins/protein_samplesheet.csv \
  --egapx_paramfile /absolute/path/to/project/egapx/input_egapx.yaml \
  --egapx_runner_dir /absolute/path/to/project/egapx_runner \
  --egapx_local_cache_dir /absolute/path/to/project/egapx_cache \
  -resume
```

Do not pass `--workflow`; TITAN rejects it because partial public modes have
been removed.

## Resume

Use `--resume` for normal restarts through the launcher. It appends Nextflow
`-resume`.

For critical production reruns, prefer pinning the exact Nextflow session UUID
when using a site launcher that supports it, for example:

```bash
TITAN_RESUME_ID=8228eeab-ca29-4701-814a-c2abc25a209c \
examples/colmar/launch_TITAN_serveur_colmar.sh
```

Plain `-resume` resumes from the latest session recorded in `.nextflow/history`.
Avoid running unrelated `nextflow run` commands from the same project directory
between production attempts, because that can make the latest session a test or
preview run.

## Outputs

Main result groups are published under `--output_dir`:

```text
aegis_outputs/
Diamond2GO_outputs/
EggNOG_outputs/
InterProScan_outputs/
additional_annotations/
final_annotations/
quality_report/
validation/
provenance/
egapx/
intermediate_files/
nextflow_reports/
```

EGAPx outputs include:

```text
egapx.complete.genomic.gff3
egapx.complete.genomic.gtf
egapx.complete.proteins.faa
egapx.complete.cds.fna
egapx.complete.transcripts.fna
egapx.annotated_genome.asn
```

`provenance/evidence_manifest.json` records main inputs, AEGIS evidence files,
final AEGIS outputs, file sizes and SHA-256 checksums. Final structural
validation reports are written under `validation/`. See
[inputs_outputs.md](inputs_outputs.md) for the complete input and output tree.

## Production Checklist

Before launching a long run:

* Confirm all TITAN and EGAPx paths are absolute and visible on compute nodes.
* Confirm the selected container runtime works on compute nodes.
* Confirm EGAPx, EDTA and AEGIS images are available or pullable.
* Pre-stage the EGAPx runner, support cache and executor config for
  reproducible HPC runs.
* Confirm the EGAPx cache contains the BUSCO lineage selected from the YAML
  `taxid`.
* Stage optional databases before enabling eggNOG-mapper, Helixer,
  InterProScan, Rfam, OMArk, CPAT or BUSCO.
* Enable tRNAscan-SE and Rfam when using the lncRNA candidate branch.
* Enable FLAIR and SQANTI3 only for runs with long-read evidence.
* Keep `run_transdecoder = true` with Mikado unless a transcript-only run is
  intentional.
* Run a `-stub-run` after every profile/config edit.
* Keep `.nextflow.log`, `nextflow_reports/`, trace/timeline files when enabled
  and the published `provenance/` directory with the run outputs.

For troubleshooting details and current limitations, see the README sections
`Troubleshooting` and `Limitations`.
