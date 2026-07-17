# TITAN installation and input preparation

This guide describes a production-oriented TITAN setup from a fresh clone to a runnable command. TITAN has one public execution contract: it always runs evidence generation, mandatory EDTA, mandatory EGAPx and Aegis in the same Nextflow graph.

## 1. System requirements

Use a Linux host or HPC login node with access to compute nodes that can run containers.

Required software:

* Git.
* Bash.
* Java compatible with Nextflow.
* Nextflow. TITAN currently declares `nextflowVersion = 24.04.3`; current newer Nextflow versions may run but can produce warnings.
* Python 3.11+ with `pyyaml`. The EGAPx runner is a Python script and needs YAML support.
* `curl` and `tar`, used by the EGAPx module to download the official EGAPx runner when `egapx_runner_dir` is not provided.
* Docker for local/containerized runs, or Apptainer/Singularity for HPC runs.

Typical checks:

```bash
git --version
java -version
nextflow -version
python3 -c 'import yaml; print("pyyaml OK")'
curl --version
tar --version
docker --version        # local Docker runs
apptainer --version     # HPC Apptainer runs
```

Install `pyyaml` if missing:

```bash
python3 -m pip install --user pyyaml
```

## 2. Clone TITAN

Clone the workflow repository and move to the TITAN directory:

```bash
git clone git@github-amandine:Grapedia/workflows.git
cd workflows/TITAN
```

If you do not have SSH access configured for the GitHub remote, use the HTTPS URL from the repository instead.

## 3. Container images

TITAN process modules use digest-pinned container images declared in `nextflow.config` as `params.container_*`. Pull at least the mandatory EGAPx and AEGIS images before a production run:

```bash
docker pull ncbi/egapx@sha256:bc657b232d93364d5f3b75ad3bfaf14b6267e46173672b609f26078d48a04298
docker pull tomsbiolab/aegis@sha256:de88470b3fb4fbab3ff2d5fa0fb9fed36b55952d1e383d3fdb2f5a3a530d84e6
```

Validate the runtime container contract after edits:

```bash
python3 scripts/validate_container_pins.py
python3 scripts/validate_profiles.py
```

For Apptainer/Singularity environments, pre-pull the same Docker image if your cluster does not allow runtime image downloads:

```bash
apptainer pull egapx_0.5.2.sif docker://ncbi/egapx@sha256:bc657b232d93364d5f3b75ad3bfaf14b6267e46173672b609f26078d48a04298
apptainer pull aegis_v0.3.25.sif docker://tomsbiolab/aegis@sha256:de88470b3fb4fbab3ff2d5fa0fb9fed36b55952d1e383d3fdb2f5a3a530d84e6
apptainer pull edta_2.2.0.sif docker://quay.io/biocontainers/edta@sha256:793cbb17bc0569e01caa0c83ad8d1756a394c2ee47b3f512ad4077bc3e422579
```

EGAPx is a nested Nextflow workflow. The TITAN EGAPx process runs the official EGAPx runner `v0.5.2`, and that runner launches EGAPx tasks using `params.egapx_executor` and `params.container_egapx`. The `apptainer` profile sets `egapx_executor = singularity` for the nested EGAPx run.

EDTA uses the official BioContainers image corresponding to tag `2.2.0--hdfd78af_1`, pinned by digest. For direct HPC testing outside TITAN:

```bash
SINGULARITY_CACHEDIR=./
export SINGULARITY_CACHEDIR
unset -f which
singularity pull EDTA.sif docker://quay.io/biocontainers/edta:2.2.0--hdfd78af_1
export PYTHONNOUSERSITE=1
singularity exec EDTA.sif EDTA.pl --genome genome.fa --species others --threads 8
```

Defaults:

```text
egapx_version = 0.5.2
egapx_revision = v0.5.2
container_egapx = ncbi/egapx@sha256:bc657b232d93364d5f3b75ad3bfaf14b6267e46173672b609f26078d48a04298
container_edta = quay.io/biocontainers/edta@sha256:793cbb17bc0569e01caa0c83ad8d1756a394c2ee47b3f512ad4077bc3e422579
egapx_executor = docker
egapx_data_version = current_1
aegis_version = v0.3.25
container_aegis = tomsbiolab/aegis@sha256:de88470b3fb4fbab3ff2d5fa0fb9fed36b55952d1e383d3fdb2f5a3a530d84e6
```

AEGIS now uses the upstream CLI container and runs `aegis merge` followed by `aegis extract`. On HPC without Docker, set `--egapx_executor singularity` for nested EGAPx execution and make the AEGIS image available through your site container runtime. Legacy `--egapx_container` and `--aegis_container` overrides are still accepted, but new runs should use `--container_egapx` and `--container_aegis`.

## 4. Required input files

Prepare these files before running TITAN:

| TITAN parameter | Required content |
| --- | --- |
| `--new_assembly` | Target genome assembly to annotate, FASTA. |
| `--previous_assembly` | Previous/reference assembly used by Liftoff, FASTA. |
| `--previous_annotations` | Annotation on the previous assembly, GFF3. |
| `--RNAseq_samplesheet` | CSV with RNA-seq inputs. |
| `--RNAseq_data_dir` | Directory containing local RNA-seq FASTQ/FASTA files referenced by the RNA-seq samplesheet. |
| `--protein_samplesheet` | CSV listing protein FASTA evidence files. |
| `--egapx_paramfile` | EGAPx YAML input file. |
| `--output_dir` | Output directory for TITAN results. |

Recommended production layout:

```text
project/
  assemblies/
    previous.fa
    target.fa
  annotations/
    previous.gff3
  rnaseq/
    RNAseq_samplesheet.csv
    sampleA.fastq.gz
    sampleB_R1.fastq.gz
    sampleB_R2.fastq.gz
    longSample.fastq.gz
  proteins/
    protein_samplesheet.csv
    araport.fa
    swissprot_plants.fa
  egapx/
    input_egapx.yaml
  titan_out/
```

## 5. RNA-seq samplesheet

TITAN expects a CSV with this header:

```csv
sampleID,SRA_or_FASTQ,library_layout
```

Allowed `library_layout` values:

* `single`: single-end short reads.
* `paired`: paired-end short reads.
* `long`: long-read RNA-seq.

Allowed `SRA_or_FASTQ` values:

* `FASTQ`: local gzip-compressed FASTQ files under `RNAseq_data_dir`.
* `FASTA`: local FASTA file under `RNAseq_data_dir`; only valid with `library_layout=long`.
* `SRA`: SRR, ERR or DRR run accession matching `sampleID`; TITAN resolves ENA FASTQ URLs, downloads the gzipped FASTQ files, retries transient failures and verifies ENA MD5 checksums when available.

Long-read processing is automatic: if at least one row has `library_layout=long`, TITAN launches the long-read branch and Aegis receives long-read evidence.

Example:

```csv
sampleID,SRA_or_FASTQ,library_layout
leaf_single,FASTQ,single
berry_paired,FASTQ,paired
isoseq_leaf,FASTQ,long
```

Local files are inferred from `sampleID` under `RNAseq_data_dir`. Keep file names consistent with existing module expectations:

```text
single: <sampleID>.fastq.gz
paired: <sampleID>_1.fastq.gz and <sampleID>_2.fastq.gz
long:   <sampleID>.fastq.gz or <sampleID>.fasta
```

SRA accessions are accepted when the run is available through ENA FASTQ metadata. Production runs are still more reproducible when FASTQ/FASTA files are downloaded and versioned outside TITAN.

## 6. Protein samplesheet

TITAN expects a CSV with this header:

```csv
organism,filename
```

Example:

```csv
organism,filename
Araport,araport.fa
SwissProtPlants,swissprot_plants.fa
```

The `filename` entries should point to FASTA protein files available to BRAKER3. Relative paths are resolved from the TITAN project directory.

Protein files must be multi-FASTA protein files suitable for ProtHint: every record needs a header, a non-empty sequence, and a unique first-word sequence ID. TITAN normalizes headers before BRAKER3 and cleans terminal `*`, `.` and `-`, but it fails early on internal `*` stop codons or non-protein characters because those usually indicate an input formatting or translation problem.

## 7. EGAPx YAML input

EGAPx takes its own YAML file through `--egapx_paramfile`. TITAN passes that file to the official EGAPx runner.

Minimum fields expected by EGAPx:

```yaml
genome: /absolute/path/to/project/assemblies/target.fa
taxid: 29760
organism: Vitis vinifera
```

TITAN validates that the YAML contains `genome`, `taxid` and `organism`, and that `genome` points to a readable FASTA file. Relative `genome` paths are resolved from the YAML file directory.

TITAN production runs should provide RNA-seq evidence to EGAPx as well:

```yaml
genome: /absolute/path/to/project/assemblies/target.fa
taxid: 29760
short_reads:
  - - leaf_single
    - - /absolute/path/to/project/rnaseq/leaf_single.fastq.gz
  - - berry_paired
    - - /absolute/path/to/project/rnaseq/berry_paired_1.fastq.gz
      - /absolute/path/to/project/rnaseq/berry_paired_2.fastq.gz
```

Recommended production example with long reads and a locus tag prefix:

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

* `genome` should normally be an absolute path visible from the node running the EGAPx runner.
* `taxid` must be the NCBI taxonomy identifier for the target organism.
* RNA-seq is highly recommended for EGAPx annotation quality.
* Local RNA-seq files should be FASTA or FASTQ, not BAM.
* If using SRA accessions, provide one run accession per row and check ENA availability before launching production.
* EGAPx performs its own masking; it does not require the EDTA-masked genome as input.
* EGAPx currently supports many vertebrates, arthropods, echinoderms and plants; fungi, protists and nematodes are out of scope according to NCBI documentation.
* TITAN's `egapx` process launches a nested EGAPx Nextflow workflow from the host environment. The host must provide `python3`, `curl`, `tar` and the selected nested executor (`docker`, `singularity` or `apptainer`).
* For strict offline reproducibility, pre-stage the EGAPx runner matching `egapx_revision` and set `--egapx_runner_dir`; otherwise TITAN downloads that pinned GitHub revision inside the task.
* Nested EGAPx work and logs live under the task work directory in `egapx_work` and under published `egapx_out/nextflow`; the outer TITAN process only stages the final named EGAPx outputs.

TITAN publishes EGAPx outputs under `${output_dir}/egapx`, including:

```text
egapx.complete.genomic.gff3
egapx.complete.genomic.gtf
egapx.complete.proteins.faa
egapx.complete.cds.fna
egapx.complete.transcripts.fna
egapx.annotated_genome.asn
egapx_out/
versions.yml
```

`egapx.complete.genomic.gff3` is passed to AEGIS as an additional annotation evidence during the final merge.

TITAN also publishes provenance under `${output_dir}/provenance`:

```text
evidence_manifest.json
versions.yml
```

`evidence_manifest.json` records the main inputs, AEGIS evidence files, final AEGIS outputs, file sizes and SHA-256 checksums.

TITAN publishes final structural validation reports under `${output_dir}/validation`:

```text
final_annotation_validation.json
final_annotation_validation.txt
```

This validation checks final GFF3 structure, coordinates, seqids, Parent links, CDS phase and protein FASTA integrity. Critical errors fail the workflow.

## 8. Validate with the built-in test profile

Before running production, verify the local workflow bootstrap:

```bash
scripts/run-tests.sh
```

This validates input schemas, parameter wiring, profile resolution, container pinning, fixture integrity and channel contracts in stub mode. It does not validate biological output quality.

## 9. Production launch

The recommended entry point is `launch_TITAN_example.sh`, which validates inputs and generates Nextflow reports:

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

The script also accepts environment variables such as `TITAN_OUTPUT_DIR`, `TITAN_NEW_ASSEMBLY`, `TITAN_EGAPX_PARAMFILE`, and `TITAN_PROFILE`.

`launch_TITAN_example.sh` resolves input paths to absolute paths, validates the container/profile contracts, and creates `TITAN_APPTAINER_CACHEDIR` automatically under the output directory when the selected profile includes `apptainer`.

Profile conventions:

```bash
nextflow config -profile local
nextflow config -profile apptainer
nextflow config -profile slurm,apptainer
nextflow config -profile test,slurm,apptainer
```

Use `test,slurm,apptainer` only for configuration resolution on the minimal fixtures. Production HPC runs should use `slurm,apptainer`.

Runtime resources are centralized in `conf/base.config` through process labels such as `process_index`, `process_alignment`, `process_transcriptome`, `process_prediction`, `process_merge` and `process_aegis`. Adjust resources in config profiles or with Nextflow process selectors rather than editing module files.

For a direct Nextflow command:

```bash
nextflow run main.nf \
  -profile slurm,apptainer \
  -name titan_target_v1 \
  -with-report titan_target_v1.report.html \
  -with-timeline titan_target_v1.timeline.html \
  -with-trace titan_target_v1.trace.txt \
  -with-dag titan_target_v1.dag.html \
  --output_dir /absolute/path/to/project/titan_out \
  --previous_assembly /absolute/path/to/project/assemblies/previous.fa \
  --new_assembly /absolute/path/to/project/assemblies/target.fa \
  --previous_annotations /absolute/path/to/project/annotations/previous.gff3 \
  --RNAseq_samplesheet /absolute/path/to/project/rnaseq/RNAseq_samplesheet.csv \
  --RNAseq_data_dir /absolute/path/to/project/rnaseq \
  --protein_samplesheet /absolute/path/to/project/proteins/protein_samplesheet.csv \
  --egapx_paramfile /absolute/path/to/project/egapx/input_egapx.yaml \
  -resume
```

Do not pass `--workflow`; TITAN rejects it because partial public modes have been removed.

## 10. Production checklist

Before launching a long run:

* Confirm all paths in TITAN params and EGAPx YAML are absolute and visible on compute nodes.
* Confirm the container runtime works on compute nodes.
* Confirm EGAPx can pull or access the pinned `ncbi/egapx@sha256:...` image from `nextflow.config`.
* Confirm EDTA can pull or access the pinned `quay.io/biocontainers/edta@sha256:...` image from `nextflow.config`.
* Confirm AEGIS can pull or access the pinned `tomsbiolab/aegis@sha256:...` image from `nextflow.config`.
* Confirm `TITAN_APPTAINER_CACHEDIR` points to a writable shared filesystem when using `slurm,apptainer`.
* Run a `-stub-run` after every config/profile edit.
* Use `-resume` for restart after interrupted runs.
* Keep `nextflow_reports/`, `trace.txt`, `timeline.html` and `.nextflow.log` with the run outputs.

For quick troubleshooting and the current limitations of stub tests, CI, SRA handling and cluster-specific Apptainer behavior, see the README sections `Troubleshooting` and `Limitations`.

## 11. References

* EGAPx official repository and `v0.5.2` runner: <https://github.com/ncbi/egapx/tree/v0.5.2>
* EGAPx input format and requirements: <https://github.com/ncbi/egapx/blob/v0.5.2/README.md>
* EGAPx official Docker image: <https://hub.docker.com/r/ncbi/egapx>
* EDTA BioContainers image: <https://quay.io/repository/biocontainers/edta>
* AEGIS repository and CLI documentation: <https://github.com/Tomsbiolab/aegis>
* AEGIS Docker image: <https://hub.docker.com/r/tomsbiolab/aegis>
