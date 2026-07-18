# TITAN installation

This page is the short setup path from a fresh clone to a runnable TITAN
command. Data formats, reference database downloads and production operations
are documented separately:

* [Input preparation](inputs.md)
* [Reference data](reference-data.md)
* [Production run guide](production-run.md)
* [Tool reference](../reference/tools.md)

TITAN has one public execution contract: evidence generation, mandatory EDTA,
mandatory EGAPx and AEGIS run in the same Nextflow graph.

## 1. Requirements

Use a Linux workstation or HPC login node with access to compute nodes that can
run containers.

Required software:

* Git.
* Bash.
* Java compatible with Nextflow.
* Nextflow. TITAN declares `nextflowVersion = 24.04.3`; newer versions may run
  but can emit warnings.
* Python 3.11+ with `pyyaml`.
* `curl` and `tar`; EGAPx runner preparation uses both.
* Docker for local runs, or Apptainer/Singularity for HPC runs.

Check the environment:

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

Install `pyyaml` if needed:

```bash
python3 -m pip install --user pyyaml
```

## 2. Clone

```bash
git clone git@github.com:Grapedia/workflows.git
cd workflows/TITAN
```

If SSH access is not configured for GitHub, use the repository HTTPS URL
instead.

## 3. Validate the Checkout

Run the bundled checks before preparing a production run:

```bash
scripts/run-tests.sh
```

This validates input schemas, parameter wiring, profile resolution, container
pinning, minimal fixture integrity and channel contracts in stub mode. It does
not validate biological output quality.

Useful individual checks after config edits:

```bash
python3 scripts/validate_container_pins.py
python3 scripts/validate_profiles.py
python3 scripts/validate_nextflow_quality.py
nextflow config -profile local
nextflow config -profile slurm,apptainer
```

## 4. Container Runtime

TITAN process modules use digest-pinned images declared in `nextflow.config` as
`params.container_*`. Docker can pull these directly. Apptainer/Singularity can
pull the same references with `docker://`.

For restricted HPC environments, pre-pull or cache at least the mandatory
images before production:

```bash
apptainer pull egapx_0.5.2.sif docker://ncbi/egapx@sha256:bc657b232d93364d5f3b75ad3bfaf14b6267e46173672b609f26078d48a04298
apptainer pull aegis_v0.3.25.sif docker://tomsbiolab/aegis@sha256:de88470b3fb4fbab3ff2d5fa0fb9fed36b55952d1e383d3fdb2f5a3a530d84e6
apptainer pull edta_2.2.0.sif docker://quay.io/biocontainers/edta@sha256:793cbb17bc0569e01caa0c83ad8d1756a394c2ee47b3f512ad4077bc3e422579
```

The complete lock list is maintained in
[container-locks.md](../development/container-locks.md).

## 5. Prepare Inputs

Required production inputs are:

| TITAN parameter | Content |
| --- | --- |
| `--new_assembly` | Target genome assembly FASTA. |
| `--previous_assembly` | Previous/reference assembly FASTA for Liftoff. |
| `--previous_annotations` | GFF3 matching `--previous_assembly`. |
| `--RNAseq_samplesheet` | CSV with RNA-seq inputs. |
| `--RNAseq_data_dir` | Directory containing local RNA-seq files. |
| `--protein_samplesheet` | CSV listing protein FASTA evidence files. |
| `--egapx_paramfile` | EGAPx YAML input file. |
| `--output_dir` | Output directory for TITAN results. |

See [inputs.md](inputs.md) for samplesheet formats, naming rules and EGAPx YAML
examples.

## 6. Prepare Reference Data

Mandatory EGAPx support data and optional datasets should be staged before long
HPC runs. The launcher can prepare EGAPx, eggNOG, Helixer, InterProScan, Rfam,
OMArk and CPAT data; BUSCO lineage downloads are prepared outside TITAN with
BUSCO's own downloader.

See [reference-data.md](reference-data.md) for the full command list and the
parameters passed back to TITAN.

## 7. First Launch

Use `launch_TITAN_example.sh` for production-oriented runs. It resolves paths,
validates inputs, checks profile/container contracts, prepares supported
reference data when requested and then launches Nextflow.

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

Use `--dry-run` to print the generated Nextflow command without executing it:

```bash
./launch_TITAN_example.sh ... --dry-run
```

## 8. Next Steps

For complete production usage, including optional branches, resume discipline,
Slurm/Apptainer conventions and final checks, continue with
[production-run.md](production-run.md).
