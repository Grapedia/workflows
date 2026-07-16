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

TITAN process modules use container images declared in `modules/*.nf`. Pull at least the mandatory EGAPx and AEGIS images before a production run:

```bash
docker pull ncbi/egapx:0.5.2
docker image inspect ncbi/egapx:0.5.2 --format '{{index .RepoDigests 0}}'
docker pull tomsbiolab/aegis:latest
docker image inspect tomsbiolab/aegis:latest --format '{{index .RepoDigests 0}}'
```

The expected image digests at the time of integration were:

```text
ncbi/egapx@sha256:bc657b232d93364d5f3b75ad3bfaf14b6267e46173672b609f26078d48a04298
tomsbiolab/aegis@sha256:de88470b3fb4fbab3ff2d5fa0fb9fed36b55952d1e383d3fdb2f5a3a530d84e6
```

For Apptainer/Singularity environments, pre-pull the same Docker image if your cluster does not allow runtime image downloads:

```bash
apptainer pull egapx_0.5.2.sif docker://ncbi/egapx:0.5.2
apptainer pull aegis_latest.sif docker://tomsbiolab/aegis:latest
```

EGAPx is a nested Nextflow workflow. The TITAN EGAPx process runs the official EGAPx runner `v0.5.2`, and that runner launches EGAPx tasks using `params.egapx_executor` and `params.egapx_container`.

Defaults:

```text
egapx_version = 0.5.2
egapx_revision = v0.5.2
egapx_container = ncbi/egapx:0.5.2
egapx_executor = docker
egapx_data_version = current_1
aegis_version = v0.3.25
aegis_container = tomsbiolab/aegis:latest
```

AEGIS now uses the upstream CLI container and runs `aegis merge` followed by `aegis extract`. On HPC without Docker, set `--egapx_executor singularity` for nested EGAPx execution and make the AEGIS image available through your site container runtime.

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

Long-read processing is automatic: if at least one row has `library_layout=long`, TITAN launches the long-read branch and Aegis receives long-read evidence.

Example:

```csv
sampleID,SRA_or_FASTQ,library_layout
leaf_single,leaf_single.fastq.gz,single
berry_paired,berry_paired,paired
isoseq_leaf,isoseq_leaf.fastq.gz,long
```

Current historical modules infer local files from `sampleID` under `RNAseq_data_dir`. Keep file names consistent with existing module expectations:

```text
single: <sampleID>.fastq.gz
paired: <sampleID>_1.fastq.gz and <sampleID>_2.fastq.gz
long:   <sampleID>.fastq.gz or <sampleID>.fasta
```

SRA accessions are still accepted by the historical preparation modules, but production runs are more reproducible when FASTQ/FASTA files are downloaded and versioned outside TITAN.

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

The `filename` entries should point to FASTA protein files available to BRAKER3. Some current modules still mount `${projectDir}/data/protein_data`; until that is fully refactored, keep a production-compatible protein data layout or test your profile on a small run first.

## 7. EGAPx YAML input

EGAPx takes its own YAML file through `--egapx_paramfile`. TITAN passes that file to the official EGAPx runner.

Minimum fields expected by EGAPx:

```yaml
genome: /absolute/path/to/project/assemblies/target.fa
taxid: 29760
```

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
* If using SRA accessions, check the number of runs before launching production. Large SRA queries can expand to many runs.
* EGAPx performs its own masking; it does not require the EDTA-masked genome as input.
* EGAPx currently supports many vertebrates, arthropods, echinoderms and plants; fungi, protists and nematodes are out of scope according to NCBI documentation.

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

## 8. Validate with the built-in test profile

Before running production, verify the local workflow bootstrap:

```bash
nextflow config -profile test
python3 scripts/validate_minimal_test_data.py
nextflow run main.nf -profile test -stub-run -ansi-log false
```

This validates parameter wiring and channel contracts. It does not validate biological output quality.

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
* Confirm EGAPx can pull or access `ncbi/egapx:0.5.2`.
* Confirm AEGIS can pull or access `tomsbiolab/aegis:latest`.
* Confirm `python3 -c 'import yaml'` works in the environment that executes the EGAPx process.
* Run a `-stub-run` after every config/profile edit.
* Use `-resume` for restart after interrupted runs.
* Keep `nextflow_reports/`, `trace.txt`, `timeline.html` and `.nextflow.log` with the run outputs.

## 11. References

* EGAPx official repository and `v0.5.2` runner: <https://github.com/ncbi/egapx/tree/v0.5.2>
* EGAPx input format and requirements: <https://github.com/ncbi/egapx/blob/v0.5.2/README.md>
* EGAPx official Docker image: <https://hub.docker.com/r/ncbi/egapx>
* AEGIS repository and CLI documentation: <https://github.com/Tomsbiolab/aegis>
* AEGIS Docker image: <https://hub.docker.com/r/tomsbiolab/aegis>
