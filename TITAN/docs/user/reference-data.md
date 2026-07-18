# TITAN reference data

This page lists external datasets and model downloads needed by TITAN modules.
Keep these directories on a persistent filesystem visible to compute nodes.

The container image pins themselves are documented in
[container-locks.md](../development/container-locks.md).

## Summary

| Branch | Default | Reference data | Preparation |
| --- | --- | --- | --- |
| EGAPx | on | EGAPx support cache, including selected BUSCO lineage | Launcher or pinned EGAPx runner |
| eggNOG-mapper | off | `eggnog.db`, `eggnog_proteins.dmnd`, `eggnog.taxa.db` | Bundled script |
| Helixer | off | Lineage `.h5` model | Bundled script |
| InterProScan | off | Member databases under `pfam/`, `cdd/`, `gene3d/`, ... | Bundled script |
| tRNAscan-SE | off | None beyond the container image | No data step |
| Infernal/Rfam | off | `Rfam.cm`, `Rfam.clanin`, `cmpress` index files | Bundled script |
| lncRNA/CPAT | off | Plant-LncPipe CPAT model files | Bundled script or validation-time fetch |
| Mikado + TransDecoder | off | None beyond container images and pipeline evidence | No data step |
| FLAIR | off | None beyond container image and long-read evidence | No data step |
| SQANTI3 | off | None beyond container image and long-read evidence | No data step |
| OMArk | off | OMAmer `omamer.h5` database | Bundled script |
| BUSCO | off | BUSCO lineage dataset | BUSCO downloader outside TITAN |

## EGAPx

EGAPx is mandatory. TITAN launches the official EGAPx runner at the pinned
revision declared by `egapx_revision`.

Recommended production parameters:

```text
egapx_revision = v0.5.2
egapx_data_version = current_1
egapx_runner_dir = /absolute/path/to/project/.egapx_runner
egapx_local_cache_dir = /absolute/path/to/project/.egapx_cache
egapx_config_dir = /absolute/path/to/project/egapx_config
```

The launcher can download the full EGAPx support cache:

```bash
./launch_TITAN_example.sh \
  --prepare-egapx-cache \
  --egapx-cache-dir /absolute/path/to/project/egapx_cache \
  --egapx-runner-dir /absolute/path/to/project/egapx_runner \
  ... # required TITAN inputs
```

For a species-specific cache, run the pinned EGAPx runner directly:

```bash
python3 ui/egapx.py input_egapx.yaml -dn -lc /absolute/path/to/project/egapx_cache -dv current_1
```

`-dn` uses the YAML `taxid` to choose EGAPx support files and the BUSCO lineage.
For `taxid: 29760` (`Vitis vinifera`), the current EGAPx runner selected
`busco_downloads/lineages/eudicots_odb10`.

Record the data version printed by EGAPx, for example
`egapxsupportdata_20251017`, and pass it back with `--egapx_data_version` for
reproducible reruns.

## eggNOG-mapper

eggNOG-mapper is optional and disabled by default:

```text
run_eggnog_mapper = false
eggnog_data_dir = false
```

Download the compatible database once:

```bash
scripts/download_eggnog_data.sh --data-dir /absolute/path/to/project/eggnog_data
```

The script fetches `eggnog.db`, `eggnog_proteins.dmnd` and `eggnog.taxa.db`
from `http://eggnog5.embl.de/download/emapperdb-5.0.2`, matching the pinned
eggNOG-mapper `2.1.15` container.

Enable it:

```bash
--run_eggnog_mapper true \
--eggnog_data_dir /absolute/path/to/project/eggnog_data
```

Launcher shortcut:

```bash
./launch_TITAN_example.sh --prepare-eggnog-data --enable-eggnog-mapper ... 
```

## Helixer

Helixer is optional and disabled by default:

```text
run_helixer = false
helixer_model = land_plant
helixer_model_dir = false
helixer_use_gpu = false
```

Fetch one lineage model with the container's own downloader:

```bash
scripts/download_helixer_model.sh \
  --model-dir /absolute/path/to/project/helixer_models \
  --container docker.io/gglyptodon/helixer-docker@sha256:e2294eb2c282c35b919933daa0d6c145b635bfac8b717dbff88e88accbde4303 \
  --lineage land_plant
```

Supported lineages are `vertebrate`, `land_plant`, `fungi` and `invertebrate`.
Helixer's container can fall back to CPU when no GPU is passed through.

Enable it:

```bash
--run_helixer true \
--helixer_model_dir /absolute/path/to/project/helixer_models \
--helixer_model land_plant
```

Launcher shortcut:

```bash
./launch_TITAN_example.sh --prepare-helixer-model --enable-helixer --helixer-lineage land_plant ...
```

Add `--enable-helixer-gpu` only when the process node can see a GPU.

## InterProScan

InterProScan is optional and disabled by default:

```text
run_interproscan = false
interproscan_data_dir = false
```

Download the member database bundle once:

```bash
scripts/download_interproscan_data.sh \
  --data-dir /absolute/path/to/project/interproscan_data
```

The script downloads the `5.78-109.0` data-only bundle from
`https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5`, verifies the upstream MD5
file and extracts member databases directly under the target directory.

Enable it:

```bash
--run_interproscan true \
--interproscan_data_dir /absolute/path/to/project/interproscan_data
```

Launcher shortcut:

```bash
./launch_TITAN_example.sh --prepare-interproscan-data --enable-interproscan ...
```

TITAN bind-mounts this directory to `/opt/interproscan/data` inside the
InterProScan container.

## tRNAscan-SE

tRNAscan-SE is optional and disabled by default. It requires no external data
directory beyond the pinned container image.

```bash
--run_trnascan true
```

## Infernal/Rfam

Rfam is optional and disabled by default:

```text
run_rfam = false
rfam_data_dir = false
```

Download and index the covariance model library:

```bash
scripts/download_rfam_data.sh \
  --data-dir /absolute/path/to/project/rfam_data \
  --container quay.io/biocontainers/infernal@sha256:05ae1ca6cc76c27180524bc38c5b1e17adf9377be5b8c644d3e8e707848d4d99
```

The script fetches `Rfam.cm.gz` and `Rfam.clanin` from
`https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT`, decompresses `Rfam.cm` and
runs `cmpress` in the pinned Infernal container.

Enable it:

```bash
--run_rfam true \
--rfam_data_dir /absolute/path/to/project/rfam_data
```

Launcher shortcut:

```bash
./launch_TITAN_example.sh --prepare-rfam-data --enable-rfam ...
```

## lncRNA/CPAT

The lncRNA branch is optional and disabled by default:

```text
run_lncrna = false
lncrna_require_cpat_model = true
cpat_model_flavour = plant_lncpipe
cpat_plant_cutoff = 0.46
```

Download the Plant-LncPipe CPAT model:

```bash
scripts/download_cpat_plant_lncpipe.sh \
  --model-dir /absolute/path/to/project/cpat_plant_lncpipe
```

The script downloads and verifies `Plant_Hexamer.tsv` and `Plant.logit.RData`
from the Plant-LncPipe model source.

Enable it:

```bash
--run_lncrna true \
--cpat_model_dir /absolute/path/to/project/cpat_plant_lncpipe
```

Launcher shortcut:

```bash
./launch_TITAN_example.sh --prepare-cpat-model --enable-lncrna ...
```

For the cleanest lncRNA candidate filtering, also enable tRNAscan-SE and Rfam.

## Mikado and TransDecoder

Mikado is optional and disabled by default. TransDecoder is enabled whenever
Mikado is enabled unless explicitly disabled.

```bash
--run_mikado true
```

Disable TransDecoder only for an intentional transcript-only Mikado run:

```bash
--run_mikado true \
--run_transdecoder false
```

No additional reference database is required.

## FLAIR

FLAIR is optional and disabled by default. It is useful only when the RNA-seq
samplesheet contains at least one `library_layout=long` row.

```bash
--run_flair true
```

No additional reference database is required.

## SQANTI3

SQANTI3 is optional and disabled by default. It is useful for long-read runs,
especially with FLAIR enabled.

```bash
--run_sqanti3 true
```

No additional reference database is required. The default
`sqanti3_libbz2_path` is a path inside the pinned SQANTI3 container:

```text
sqanti3_libbz2_path = /usr/local/lib/libbz2.so.1.0.8
```

If the SQANTI3 container is changed, set this parameter to the matching
container-visible library path or to `false` if the image already provides
`libbz2.so.1`.

## OMArk

OMArk is optional and disabled by default:

```text
run_omark = false
omark_data_dir = false
```

Download the OMAmer database:

```bash
scripts/download_omark_data.sh --data-dir /absolute/path/to/project/omark_data
```

The default URL is `https://omabrowser.org/All/LUCA.h5`. The file is large and
is saved as `omamer.h5`.

Enable it:

```bash
--run_omark true \
--omark_data_dir /absolute/path/to/project/omark_data
```

Launcher shortcut:

```bash
./launch_TITAN_example.sh --prepare-omark-data --enable-omark ...
```

## BUSCO

BUSCO is optional and disabled by default:

```text
run_busco = false
busco_lineage = eudicotyledons_odb12.2
busco_data_dir = false
```

TITAN does not bundle a BUSCO downloader. Prepare the selected lineage outside
TITAN, for example:

```bash
busco --download eudicotyledons_odb12.2 --download_path /absolute/path/to/project/busco_data
```

Enable it:

```bash
--run_busco true \
--busco_data_dir /absolute/path/to/project/busco_data \
--busco_lineage eudicotyledons_odb12.2
```

TITAN validates `--busco_data_dir` whenever `--run_busco true` is set.

## Colmar Launcher

The Colmar example launcher exposes preparation flags for the datasets TITAN
can stage automatically:

```bash
examples/colmar/launch_TITAN_serveur_colmar.sh --prepare-egapx-cache
examples/colmar/launch_TITAN_serveur_colmar.sh --prepare-eggnog-data
examples/colmar/launch_TITAN_serveur_colmar.sh --prepare-helixer-model
examples/colmar/launch_TITAN_serveur_colmar.sh --prepare-interproscan-data
examples/colmar/launch_TITAN_serveur_colmar.sh --prepare-rfam-data
examples/colmar/launch_TITAN_serveur_colmar.sh --prepare-omark-data
examples/colmar/launch_TITAN_serveur_colmar.sh --prepare-cpat-model
```

Enable optional branches in the site copy of
`examples/colmar/slurm_apptainer.config`, or pass the equivalent Nextflow
parameters directly.

## References

* EGAPx repository: <https://github.com/ncbi/egapx/tree/v0.5.2>
* eggNOG-mapper data: <http://eggnog5.embl.de/download/emapperdb-5.0.2>
* Helixer: <https://github.com/weberlab-hhu/Helixer>
* InterProScan data: <https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5>
* Rfam data: <https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT>
* Plant-LncPipe CPAT model: <https://github.com/xuechantian/Plant-LncRNA-pipline>
* OMAmer databases: <https://omabrowser.org/All/>
* BUSCO: <https://gitlab.com/ezlab/busco>
