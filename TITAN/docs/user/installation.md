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
git clone git@github.com:Grapedia/workflows.git
cd workflows/TITAN
```

If you do not have SSH access configured for the GitHub remote, use the HTTPS URL from the repository instead.

## 3. Container images

TITAN process modules use digest-pinned container images declared in `nextflow.config` as `params.container_*`. Pull at least the mandatory EGAPx and AEGIS images before a production run:

```bash
docker pull ncbi/egapx@sha256:bc657b232d93364d5f3b75ad3bfaf14b6267e46173672b609f26078d48a04298
docker pull tomsbiolab/aegis@sha256:de88470b3fb4fbab3ff2d5fa0fb9fed36b55952d1e383d3fdb2f5a3a530d84e6
```

eggNOG-mapper is optional (see [section 8](#8-eggnog-mapper-optional)); pull it only if you plan to enable `--run_eggnog_mapper true`:

```bash
docker pull quay.io/biocontainers/eggnog-mapper@sha256:f70babaf681ff4b6b2fc8e8e76754bf989f01dd4910ca91f156c22aa88ea70d3
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

AEGIS now uses the upstream CLI container and runs `aegis merge`, then `aegis rename` (gene/transcript/CDS/exon IDs regenerated under the `--aegis_gene_id_prefix` scheme, default `Vitvi`, with an old-to-new correspondences TSV), then `aegis tidy` (`--standard-features`, e.g. normalizing `pseudotranscript` to `mRNA`), and finally `aegis extract`. On HPC without Docker, set `--egapx_executor singularity` for nested EGAPx execution and make the AEGIS image available through your site container runtime. Legacy `--egapx_container` and `--aegis_container` overrides are still accepted, but new runs should use `--container_egapx` and `--container_aegis`.

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

### EGAPx support cache and BUSCO lineage

EGAPx can run online, but production HPC runs are more reproducible when its support data are pre-staged in a local cache and passed to TITAN with `--egapx_local_cache_dir`.

TITAN does not hardcode EGAPx reference-data URLs. Reference discovery and downloads are handled by the pinned EGAPx runner and the selected `--egapx_data_version`; TITAN only passes the cache directory and version through to that runner.

EGAPx has two cache download modes:

```bash
# Download the complete EGAPx support manifest for the selected data version.
python3 ui/egapx.py input_egapx.yaml -dl -lc /path/to/egapx_cache -dv current_1

# Download support data needed by this input YAML, including SRA inputs if any
# are listed as accessions, and the BUSCO lineage chosen from taxid.
python3 ui/egapx.py input_egapx.yaml -dn -lc /path/to/egapx_cache -dv current_1
```

`-dn` is preferred when preparing a cache for a single production species. It uses the YAML `taxid` to choose the closest protein set, HMM parameters, orthology reference and BUSCO lineage. It is not limited to one literal `<taxid>.params` file: EGAPx also needs common support files such as taxonomy, scoring rules, reference sets and the selected BUSCO lineage.

For `taxid: 29760` (`Vitis vinifera`), the current EGAPx runner selected the BUSCO lineage `eudicots_odb10`. A valid cache therefore contains a directory like:

```text
egapx_cache/
  busco_downloads/lineages/eudicots_odb10/
  gnomon/
  misc/
  ortholog_references/
  reference_sets/
  target_proteins/
  taxonomy/
```

If the lineage is missing, EGAPx fails before launching its nested workflow and prints a command similar to:

```bash
python3 ui/egapx.py input_egapx.yaml -lc /path/to/egapx_cache -dl
```

For smaller species-specific caches, use `-dn` instead of the suggested broad `-dl` command:

```bash
python3 ui/egapx.py input_egapx.yaml -lc /path/to/egapx_cache -dn -dv egapxsupportdata_20251017
```

Record the data version printed by EGAPx, for example `egapxsupportdata_20251017`, and pass it back to TITAN with `--egapx_data_version` for reproducible reruns.

### EGAPx executor config

EGAPx creates an `egapx_config/` directory the first time it is run with a new executor and asks the user to edit it. TITAN can pass a prepared directory with `--egapx_config_dir`.

A typical Slurm/Apptainer setup uses a local directory like:

```text
egapx_config/
  singularity.config
  process_resources.config
```

That config can run EGAPx's nested Nextflow locally inside the outer TITAN Slurm job while using Apptainer/Singularity containers. This avoids nested Slurm submission for small or moderate runs. For full production, either:

* keep this single-node nested EGAPx model and allocate enough CPUs, memory and walltime to the outer TITAN `egapx` process; or
* create a site-specific EGAPx executor config from the official EGAPx examples if you want EGAPx itself to submit nested cluster jobs.

On this cluster, use `apptainer/1.4.0-rc.2` for production validation; `apptainer/1.4.0` showed a session-directory failure on compute nodes during testing.

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

## 8. eggNOG-mapper (optional)

eggNOG-mapper adds orthology-based functional annotation on the AEGIS-derived protein FASTAs, alongside Diamond2GO. It is disabled by default (`run_eggnog_mapper = false`) because it needs an external, pre-downloaded eggNOG database that TITAN does not fetch on its own during a run.

Download the database once with the bundled script:

```bash
scripts/download_eggnog_data.sh --data-dir /absolute/path/to/project/eggnog_data
```

This fetches `eggnog.db`, `eggnog_proteins.dmnd` and `eggnog.taxa.db` with plain `curl`/`gunzip`/`tar` from `http://eggnog5.embl.de/download/emapperdb-5.0.2` (data version compatible with eggNOG-mapper `2.1.15`, matching `container_eggnog_mapper`). No eggNOG-mapper installation or container runtime is required for this step. The `download_eggnog_data.py` script bundled inside the `eggnog-mapper:2.1.15` BioContainers image points at the retired host `eggnogdb.embl.de`, which no longer resolves, so this script fetches the same data version directly from the current host instead.

Both launcher scripts can run this step automatically:

```bash
./launch_TITAN_example.sh --prepare-eggnog-data --enable-eggnog-mapper ... # plus the usual required options
examples/colmar/launch_TITAN_serveur_colmar.sh --prepare-eggnog-data
```

`--prepare-eggnog-data` downloads the database into `TITAN_EGGNOG_DATA_DIR` (default `<project-dir>/.eggnog_data`) if it is not already there; `--enable-eggnog-mapper` (in `launch_TITAN_example.sh`) then passes `--run_eggnog_mapper true --eggnog_data_dir "$TITAN_EGGNOG_DATA_DIR"` to Nextflow. On the Colmar launcher, set `run_eggnog_mapper = true` in the site copy of `examples/colmar/slurm_apptainer.config` once the database has been downloaded.

Defaults:

```text
run_eggnog_mapper = false
eggnog_data_dir = false
eggnog_mapper_cpus = 5
eggnog_mapper_sensmode = sensitive
eggnog_mapper_tax_scope = false
container_eggnog_mapper = quay.io/biocontainers/eggnog-mapper@sha256:f70babaf681ff4b6b2fc8e8e76754bf989f01dd4910ca91f156c22aa88ea70d3
```

TITAN validates `--eggnog_data_dir` before heavy execution whenever `--run_eggnog_mapper true` is set, and fails fast if the directory is missing. Outputs are published under `${output_dir}/EggNOG_outputs`.

## 9. Helixer (optional)

Helixer is an ab initio/deep-learning gene predictor run directly on the EDTA soft-masked genome (`edta.MAKER.masked`). Its GFF3 is published separately under `additional_annotations/`, and is also passed to AEGIS as optional evidence: when `run_helixer = true` it is merged into `final_annotation.gff3` alongside the other named evidence channels (renamed under the same `--aegis_gene_id_prefix` scheme as the rest of the merged annotation). It is disabled by default (`run_helixer = false`).

Helixer's official container only ships CUDA-bundled images; there is no separate CPU-only tag. This is not a problem: TensorFlow falls back to CPU automatically whenever no GPU is passed through to the container, so the same image is used for both CPU and GPU runs. GPU/CPU selection is a container-runtime flag, not a `Helixer.py` option — TITAN sets `containerOptions` per profile (`--nv` for Apptainer/Singularity, `--gpus all` for Docker) only when `--helixer_use_gpu true` is set, and only actually helps if a GPU is visible on the node.

Helixer's own `Helixer.py` does **not** auto-download its lineage model and fails fast with a clear error if it is missing (`Cannot continue without a model, either download models with fetch_helixer_models.py or set --model-filepath`). Fetch a model once with the bundled script, which runs the container's own `fetch_helixer_models.py`:

```bash
scripts/download_helixer_model.sh \
  --model-dir /absolute/path/to/project/helixer_models \
  --container docker.io/gglyptodon/helixer-docker@sha256:e2294eb2c282c35b919933daa0d6c145b635bfac8b717dbff88e88accbde4303 \
  --lineage land_plant
```

This populates `<model-dir>/<lineage>/*.h5`. The script skips the fetch if a model is already present, so re-running it is a no-op. Both launcher scripts can run this step automatically:

```bash
./launch_TITAN_example.sh --prepare-helixer-model --enable-helixer ... # plus the usual required options
examples/colmar/launch_TITAN_serveur_colmar.sh --prepare-helixer-model
```

`--prepare-helixer-model` downloads the model into `TITAN_HELIXER_MODEL_DIR` (default `<project-dir>/.helixer_models`) if it is not already there; `--enable-helixer` (in `launch_TITAN_example.sh`) then passes `--run_helixer true --helixer_model_dir "$TITAN_HELIXER_MODEL_DIR" --helixer_model "$TITAN_HELIXER_LINEAGE"` to Nextflow, and `--enable-helixer-gpu` additionally passes `--helixer_use_gpu true`. On the Colmar launcher, set `run_helixer = true` and optionally `helixer_use_gpu = true` in the site copy of `examples/colmar/slurm_apptainer.config` once the model has been downloaded.

Defaults:

```text
run_helixer = false
helixer_model_dir = false
helixer_model = land_plant
helixer_use_gpu = false
helixer_cpus = 5
container_helixer = docker.io/gglyptodon/helixer-docker@sha256:e2294eb2c282c35b919933daa0d6c145b635bfac8b717dbff88e88accbde4303
```

TITAN validates `--helixer_model_dir` and the presence of the requested `--helixer_model` lineage subdirectory before heavy execution whenever `--run_helixer true` is set, and fails fast otherwise. Outputs are published under `${output_dir}/additional_annotations/helixer`, and a dedicated `additional_annotations_manifest.json` is written under `${output_dir}/provenance` alongside the main `evidence_manifest.json`.

## 10. InterProScan (optional)

InterProScan adds member-database functional annotation (Pfam, CDD, Gene3D, PANTHER, SMART, ProSiteProfiles, ...) on the AEGIS-derived protein FASTAs, alongside Diamond2GO and eggNOG-mapper. It is disabled by default (`run_interproscan = false`) because it needs an external, pre-downloaded member database data directory that TITAN does not fetch on its own during a run.

Download the data once with the bundled script:

```bash
scripts/download_interproscan_data.sh --data-dir /absolute/path/to/project/interproscan_data
```

This fetches the official data-only bundle (`interproscan-data-5.78-109.0.tar.gz`, ~7 GB compressed) from `https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.78-109.0/alt/`, verifies it against the upstream `.md5` checksum, and extracts it so `<data-dir>` directly holds `pfam/`, `cdd/`, `gene3d/`, ... (no container or InterProScan installation required for this step). The version must match `container_interproscan`. The script skips the download if `<data-dir>/pfam` already looks populated, so re-running it is a no-op once the data is fetched.

Both launcher scripts can run this step automatically:

```bash
./launch_TITAN_example.sh --prepare-interproscan-data --enable-interproscan ... # plus the usual required options
examples/colmar/launch_TITAN_serveur_colmar.sh --prepare-interproscan-data
```

`--prepare-interproscan-data` downloads the data into `TITAN_INTERPROSCAN_DATA_DIR` (default `<project-dir>/.interproscan_data`) if it is not already there; `--enable-interproscan` (in `launch_TITAN_example.sh`) then passes `--run_interproscan true --interproscan_data_dir "$TITAN_INTERPROSCAN_DATA_DIR"` to Nextflow. On the Colmar launcher, set `run_interproscan = true` in the site copy of `examples/colmar/slurm_apptainer.config` once the data has been downloaded.

InterProScan hardcodes `data/` as a path relative to its install directory rather than accepting a CLI flag, so TITAN bind-mounts `interproscan_data_dir` to the fixed container path `/opt/interproscan/data` (matching the official Docker usage pattern) via `containerOptions` per profile, the same mechanism used for `eggnog_data_dir` and `helixer_model_dir`.

By default TITAN runs every freely available analysis (no `-appl` restriction) with `-goterms -pathways` for GO term/pathway lookups and `-dp` to disable the online precalculated match lookup service, so runs stay fully offline. This is the heaviest of the three functional annotation steps; consider it in resource planning (`process_aegis` label resources, or a dedicated `withName: interproscan` override in your production config) for large protein sets.

Defaults:

```text
run_interproscan = false
interproscan_data_dir = false
interproscan_cpus = 5
container_interproscan = docker.io/interpro/interproscan@sha256:dc58b7c147fbbf00c2dd4f5ced42121fc1e8841fcbc7cc2c484380248ff76d11
```

TITAN validates `--interproscan_data_dir` before heavy execution whenever `--run_interproscan true` is set (directory exists and contains a `pfam` subdirectory), and fails fast otherwise. Outputs are published under `${output_dir}/InterProScan_outputs`.

## 10a. tRNAscan-SE (optional)

tRNAscan-SE adds tRNA annotation directly from `--new_assembly`. It is disabled by default (`run_trnascan = false`) and does not require an offline database beyond the pinned container image.

Enable it in a direct Nextflow command with:

```bash
--run_trnascan true
```

On the Colmar production launcher, set `run_trnascan = true` in the site copy of `examples/colmar/slurm_apptainer.config`. No `--prepare-*` step is needed.

Defaults:

```text
run_trnascan = false
container_trnascan = quay.io/biocontainers/trnascan-se@sha256:e573090368974ff1228e6894828c6c8a132dfecc3198f5e9fb76832f8f434f29
```

Outputs are published under `${output_dir}/additional_annotations/ncrna/trna`. The resulting GFF3 is used by the lncRNA candidate filter and ncRNA quality summaries, but TITAN does not merge tRNAs into the final AEGIS coding annotation. See the [tRNAscan-SE tool reference](../reference/tools.md#trnascan-se).

## 10b. Infernal/Rfam ncRNA (optional)

Infernal/Rfam adds ncRNA family annotation directly from `--new_assembly`. It is disabled by default (`run_rfam = false`) because it needs a local Rfam covariance-model library.

Download and index the Rfam data once with the bundled script:

```bash
scripts/download_rfam_data.sh \
  --data-dir /absolute/path/to/project/rfam_data \
  --container quay.io/biocontainers/infernal@sha256:05ae1ca6cc76c27180524bc38c5b1e17adf9377be5b8c644d3e8e707848d4d99
```

The script downloads `Rfam.cm.gz` and `Rfam.clanin` from the current Rfam FTP directory, decompresses `Rfam.cm`, and runs `cmpress` inside the pinned Infernal container. It skips work when the indexed files are already present.

Enable it in a direct Nextflow command with:

```bash
--run_rfam true \
--rfam_data_dir /absolute/path/to/project/rfam_data
```

The Colmar production launcher can run the preparation step:

```bash
examples/colmar/launch_TITAN_serveur_colmar.sh --prepare-rfam-data
```

For that launcher, keep `run_rfam = true` and `rfam_data_dir = "${projectDir}/.rfam_data"` in the site copy of `examples/colmar/slurm_apptainer.config`, or override `TITAN_RFAM_DATA_DIR` before launching.

Defaults:

```text
run_rfam = false
rfam_data_dir = false
container_infernal = quay.io/biocontainers/infernal@sha256:05ae1ca6cc76c27180524bc38c5b1e17adf9377be5b8c644d3e8e707848d4d99
```

TITAN validates that `--rfam_data_dir` contains `Rfam.cm` and `Rfam.clanin` whenever `--run_rfam true` is set. Outputs are published under `${output_dir}/additional_annotations/ncrna/rfam`, used by lncRNA filtering and ncRNA quality summaries, and kept separate from the final AEGIS coding annotation. See the [Infernal/Rfam tool reference](../reference/tools.md#infernalrfam-ncrna).

## 10c. lncRNA/CPAT candidates (optional)

The lncRNA branch builds candidate non-coding transcripts from merged transcript evidence after AEGIS, tRNAscan-SE and Infernal/Rfam have completed. It is disabled by default (`run_lncrna = false`) and is best used together with `--run_trnascan true` and `--run_rfam true`, because those annotations help remove tRNA and known ncRNA overlaps from the candidate set.

TITAN uses the Plant-LncPipe CPAT-plant model by default. The default model directory is tracked as:

```text
cpat_model_dir = ${projectDir}/resources/cpat_plant_lncpipe
```

If the model files are missing and `lncrna_require_cpat_model = true`, TITAN downloads them during input validation with:

```bash
scripts/download_cpat_plant_lncpipe.sh --model-dir /absolute/path/to/project/resources/cpat_plant_lncpipe
```

You can also run that command manually before production, which is preferable on restricted HPC login nodes.

Enable the branch with:

```bash
--run_lncrna true \
--cpat_model_dir /absolute/path/to/project/resources/cpat_plant_lncpipe
```

Defaults:

```text
run_lncrna = false
lncrna_require_cpat_model = true
cpat_model_dir = ${projectDir}/resources/cpat_plant_lncpipe
cpat_model_flavour = plant_lncpipe
cpat_plant_cutoff = 0.46
container_cpat = quay.io/biocontainers/cpat@sha256:87366fff67d441f64e0ac4681ccbaf1147f2c0601f3df86bb99f228d7f9a9000
```

Outputs are published under `${output_dir}/additional_annotations/ncrna/lncrna`. They are candidate lncRNA annotations and are not treated as a validated final lncRNA annotation layer. See the [lncRNA candidate tool reference](../reference/tools.md#lncrna-candidates).

## 10d. Mikado and TransDecoder (optional)

Mikado produces an alternative final annotation source from the same major evidence families as AEGIS: Liftoff, EGAPx, BRAKER3, STAR/StringTie, HISAT2/StringTie, STAR/PsiCLASS, optional long-read StringTie, optional FLAIR and optional Helixer. It is disabled by default (`run_mikado = false`).

Enable Mikado with:

```bash
--run_mikado true
```

TransDecoder is enabled by default when Mikado is enabled. Disable it only if you intentionally want a Mikado run without ORF prediction:

```bash
--run_mikado true \
--run_transdecoder false
```

Defaults:

```text
run_mikado = false
run_transdecoder = true
container_mikado = quay.io/biocontainers/mikado@sha256:dd6f5a2a2d7fdbab73c835cd0f49bd1444ecaddf8e4cd96fbf0fe24f5ecf5f22
container_transdecoder = quay.io/biocontainers/transdecoder@sha256:c70f3a30cc8f3aecccb1d8978b9a49865d3994ebd7885361ab4c9dd820bd17f5
```

Outputs are published under `${output_dir}/final_annotations/mikado`. AEGIS remains the primary final annotation path; Mikado is kept as a separate final annotation source and compared against AEGIS in the quality report when enabled. See the [Mikado tool reference](../reference/tools.md#mikado-final-annotation-source).

## 10e. FLAIR long-read isoforms (optional)

FLAIR improves long-read isoform evidence. It is disabled by default (`run_flair = false`) and only produces biological outputs when the RNA-seq samplesheet contains at least one `library_layout=long` row.

Enable it with:

```bash
--run_flair true
```

FLAIR uses Minimap2 long-read alignments and Liftoff-derived splice-junction correction evidence. Its merged isoforms are passed to AEGIS and Mikado as optional transcript evidence.

Defaults:

```text
run_flair = false
container_flair = quay.io/biocontainers/flair@sha256:187e2e22535d73ecc724afc7e474d9908b6e43a55f8588e8566db3bea2eba79e
```

Outputs are published under `${output_dir}/additional_annotations/flair`. See the [FLAIR tool reference](../reference/tools.md#flair-long-read-isoforms).

## 10f. SQANTI3 long-read isoform QC (optional)

SQANTI3 evaluates long-read transcript models against the final AEGIS annotation and target genome. It is disabled by default (`run_sqanti3 = false`) and is useful only for runs with long-read evidence and, ideally, FLAIR enabled.

Enable it with:

```bash
--run_sqanti3 true
```

When long reads are absent, or when a long-read source has no isoforms, TITAN writes zero-count sentinel summaries so downstream MultiQC aggregation remains stable.

Defaults:

```text
run_sqanti3 = false
sqanti3_libbz2_path = /usr/local/lib/libbz2.so.1.0.8
container_sqanti3 = quay.io/biocontainers/sqanti3@sha256:3bd6ec96b3f1c9cae69cfef54ba0522b7d99efa7ebb0ff6a611841aa6784f74c
```

The default `sqanti3_libbz2_path` is a path inside the pinned SQANTI3 container, not a host path. TITAN uses it only to create a task-local `libbz2.so.1` compatibility symlink for `gtfToGenePred`. If you change `container_sqanti3`, either set `--sqanti3_libbz2_path` to the matching container-visible `libbz2.so.1.*` file or set it to `false` when the image already provides `libbz2.so.1`.

Outputs are published under `${output_dir}/additional_annotations/sqanti3` and `${output_dir}/quality_report/sqanti3`. See the [SQANTI3 tool reference](../reference/tools.md#sqanti3-long-read-isoform-qc).

## 10g. OMArk (optional)

OMArk adds proteome consistency, completeness and contamination checks to the quality report. It is disabled by default (`run_omark = false`) because it requires an offline OMAmer database.

Download the OMAmer database once with:

```bash
scripts/download_omark_data.sh --data-dir /absolute/path/to/project/omark_data
```

The default script downloads the full LUCA database to `omamer.h5`. This is a large multi-GB download, but it only needs to be staged once.

Enable OMArk with:

```bash
--run_omark true \
--omark_data_dir /absolute/path/to/project/omark_data
```

The Colmar production launcher can run the preparation step:

```bash
examples/colmar/launch_TITAN_serveur_colmar.sh --prepare-omark-data
```

For that launcher, keep `run_omark = true` and `omark_data_dir = "${projectDir}/.omark_data"` in the site copy of `examples/colmar/slurm_apptainer.config`, or override `TITAN_OMARK_DATA_DIR` before launching.

Defaults:

```text
run_omark = false
omark_data_dir = false
container_omark = quay.io/biocontainers/omark@sha256:84413cc19053c5d6452fbff245c9e6980b3f16aabdf991f9e51d7b9f2e0e0843
```

TITAN validates that `--omark_data_dir` contains `omamer.h5` whenever `--run_omark true` is set. Outputs are published under `${output_dir}/quality_report/omark`. See the [quality-report tool reference](../reference/tools.md#quality-report).

## 10h. BUSCO (optional)

BUSCO adds protein-mode completeness metrics to the quality report. It is disabled by default (`run_busco = false`) and needs a local BUSCO lineage dataset.

TITAN does not provide a `--prepare-busco-data` launcher flag or bundled BUSCO data downloader. Prepare the lineage outside TITAN using BUSCO's own downloader, for example:

```bash
busco --download eudicotyledons_odb12.2 --download_path /absolute/path/to/project/busco_data
```

Then enable BUSCO with:

```bash
--run_busco true \
--busco_data_dir /absolute/path/to/project/busco_data \
--busco_lineage eudicotyledons_odb12.2
```

Defaults:

```text
run_busco = false
busco_lineage = eudicotyledons_odb12.2
busco_data_dir = false
container_busco = quay.io/biocontainers/busco@sha256:d55ad622a5cafcd63c42fc309108688ab255bb9586ee756a5149e249d418c8bd
```

TITAN requires `--busco_data_dir` whenever `--run_busco true` is set. Outputs are published under `${output_dir}/quality_report/busco`. See the [quality-report tool reference](../reference/tools.md#quality-report).

## 11. Validate with the built-in test profile

Before running production, verify the local workflow bootstrap:

```bash
scripts/run-tests.sh
```

This validates input schemas, parameter wiring, profile resolution, container pinning, fixture integrity and channel contracts in stub mode. It does not validate biological output quality.

## 12. Production launch

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

## 13. Production checklist

Before launching a long run:

* Confirm all paths in TITAN params and EGAPx YAML are absolute and visible on compute nodes.
* Confirm the container runtime works on compute nodes.
* Confirm EGAPx can pull or access the pinned `ncbi/egapx@sha256:...` image from `nextflow.config`.
* Pre-stage the EGAPx runner and support cache for reproducibility, then pass `--egapx_runner_dir`, `--egapx_local_cache_dir`, `--egapx_config_dir` and a fixed `--egapx_data_version`.
* Confirm the EGAPx cache contains the BUSCO lineage selected from the YAML `taxid` (for `taxid: 29760`, `busco_downloads/lineages/eudicots_odb10`).
* Confirm EDTA can pull or access the pinned `quay.io/biocontainers/edta@sha256:...` image from `nextflow.config`.
* Confirm AEGIS can pull or access the pinned `tomsbiolab/aegis@sha256:...` image from `nextflow.config`.
* If enabling eggNOG-mapper, run `--prepare-eggnog-data` (or `scripts/download_eggnog_data.sh`) at least once and confirm `eggnog_data_dir` points to a populated database directory.
* If enabling Helixer, run `--prepare-helixer-model` (or `scripts/download_helixer_model.sh`) at least once and confirm `helixer_model_dir` contains the requested lineage; only set `helixer_use_gpu = true` if a GPU is actually visible on the node running that process.
* If enabling InterProScan, run `--prepare-interproscan-data` (or `scripts/download_interproscan_data.sh`) at least once and confirm `interproscan_data_dir` contains a populated `pfam` subdirectory; plan for the extra runtime (all analyses, on both protein FASTAs, is the heaviest of the three functional annotation steps).
* If enabling tRNAscan-SE, confirm the pinned `trnascan-se` image is available; no offline data directory is required.
* If enabling Infernal/Rfam, run `--prepare-rfam-data` on the launcher (or `scripts/download_rfam_data.sh`) at least once and confirm `rfam_data_dir` contains `Rfam.cm`, `Rfam.clanin` and the `cmpress` index files.
* If enabling lncRNA candidates, stage the CPAT Plant-LncPipe model with `scripts/download_cpat_plant_lncpipe.sh` or confirm TITAN can download it during validation; enable tRNAscan-SE and Rfam for the cleanest candidate filtering.
* If enabling Mikado, keep `run_transdecoder = true` unless you intentionally want a transcript-only Mikado run, and plan for a separate final annotation output under `final_annotations/mikado`.
* If enabling FLAIR or SQANTI3, confirm the RNA-seq samplesheet contains at least one `library_layout=long` sample.
* If enabling OMArk, run `--prepare-omark-data` on the launcher (or `scripts/download_omark_data.sh`) at least once and confirm `omark_data_dir` contains `omamer.h5`.
* If enabling BUSCO, pre-stage the selected lineage dataset and confirm `busco_data_dir` points to the directory containing that lineage.
* Confirm `TITAN_APPTAINER_CACHEDIR` points to a writable shared filesystem when using `slurm,apptainer`.
* Run a `-stub-run` after every config/profile edit.
* Use `-resume` for restart after interrupted runs.
* Keep `nextflow_reports/`, `trace.txt`, `timeline.html` and `.nextflow.log` with the run outputs.

For quick troubleshooting and the current limitations of stub tests, CI, SRA handling and cluster-specific Apptainer behavior, see the README sections `Troubleshooting` and `Limitations`.

## 14. References

* EGAPx official repository and `v0.5.2` runner: <https://github.com/ncbi/egapx/tree/v0.5.2>
* EGAPx input format and requirements: <https://github.com/ncbi/egapx/blob/v0.5.2/README.md>
* EGAPx official Docker image: <https://hub.docker.com/r/ncbi/egapx>
* EDTA BioContainers image: <https://quay.io/repository/biocontainers/edta>
* AEGIS repository and CLI documentation: <https://github.com/Tomsbiolab/aegis>
* AEGIS Docker image: <https://hub.docker.com/r/tomsbiolab/aegis>
* eggNOG-mapper repository and usage: <https://github.com/eggnogdb/eggnog-mapper>
* eggNOG-mapper BioContainers image: <https://quay.io/repository/biocontainers/eggnog-mapper>
* InterProScan repository and usage: <https://github.com/ebi-pf-team/interproscan>
* InterProScan Docker image: <https://hub.docker.com/r/interpro/interproscan>
* InterProScan member database data download: <https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/>
* Helixer repository and usage: <https://github.com/weberlab-hhu/Helixer>
* Helixer Docker image: <https://hub.docker.com/r/gglyptodon/helixer-docker>
* tRNAscan-SE repository and usage: <https://github.com/UCSC-LoweLab/tRNAscan-SE>
* Infernal repository and usage: <http://eddylab.org/infernal/>
* Rfam database downloads: <https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/>
* CPAT repository and usage: <https://github.com/liguowang/cpat>
* Plant-LncPipe CPAT plant model source: <https://github.com/xuechantian/Plant-LncRNA-pipline>
* Mikado repository and usage: <https://github.com/EI-CoreBioinformatics/mikado>
* TransDecoder repository and usage: <https://github.com/TransDecoder/TransDecoder>
* FLAIR repository and usage: <https://github.com/BrooksLabUCSC/flair>
* SQANTI3 repository and usage: <https://github.com/ConesaLab/SQANTI3>
* OMArk repository and usage: <https://github.com/DessimozLab/OMArk>
* OMAmer database downloads: <https://omabrowser.org/All/>
* BUSCO repository and usage: <https://gitlab.com/ezlab/busco>
* MultiQC repository and usage: <https://github.com/MultiQC/MultiQC>
