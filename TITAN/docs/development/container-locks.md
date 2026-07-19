# TITAN container locks

## Table Of Contents

- [Runtime images](#runtime-images)
- [Dockerfile base images](#dockerfile-base-images)
- [Test scope](#test-scope)


TITAN runtime containers are locked in `nextflow.config` with `image@sha256` references. Process modules must use `params.container_*` values instead of hardcoded image names.

Validate the contract with:

```bash
python3 scripts/validate_container_pins.py
```

This checks:

* every required `params.container_*` value exists and uses `@sha256`;
* every `modules/*.nf` `container` directive uses a central parameter;
* every `dockerfiles/**/Dockerfile` `FROM` directive uses `@sha256`;
* active runtime files do not reference `:latest`.

## Runtime images

| Tool | Parameter | Locked image |
| --- | --- | --- |
| AEGIS | `container_aegis` | `tomsbiolab/aegis@sha256:de88470b3fb4fbab3ff2d5fa0fb9fed36b55952d1e383d3fdb2f5a3a530d84e6` |
| AGAT | `container_agat` | `quay.io/biocontainers/agat@sha256:7ea8fa5a8428758cd87e3a5dcfaf277febdfcae95cd1fe473770abf8b928ec99` |
| BRAKER3 | `container_braker3` | `avelt/braker3@sha256:e69a9aaaafa81e4da5b2bbb98ae120d873018ae40453630f60051ecd5f622c44` |
| Diamond2GO | `container_diamond2go` | `avelt/diamond2go@sha256:40f1063307f98a2357d60b306bd7d79b6088591c1613e6552613da24002e8360` |
| eggNOG-mapper | `container_eggnog_mapper` | `quay.io/biocontainers/eggnog-mapper@sha256:f70babaf681ff4b6b2fc8e8e76754bf989f01dd4910ca91f156c22aa88ea70d3` |
| EDTA | `container_edta` | `quay.io/biocontainers/edta@sha256:793cbb17bc0569e01caa0c83ad8d1756a394c2ee47b3f512ad4077bc3e422579` |
| EGAPx | `container_egapx` | `ncbi/egapx@sha256:bc657b232d93364d5f3b75ad3bfaf14b6267e46173672b609f26078d48a04298` |
| fastp | `container_fastp` | `quay.io/biocontainers/fastp@sha256:0bdf8d8254fc86dd9038551d68dbcb72562e65560b9ce0ea08c1329d2f8587b4` |
| GFFCompare | `container_gffcompare` | `avelt/gffcompare@sha256:bd411c13352a2545641c8c34b701030b3977056696607e85b4e86c876d10a82c` |
| Helixer | `container_helixer` | `docker.io/gglyptodon/helixer-docker@sha256:e2294eb2c282c35b919933daa0d6c145b635bfac8b717dbff88e88accbde4303` |
| HISAT2 | `container_hisat2` | `avelt/hisat2@sha256:022933fd0d30fe9fdfd83c175f7e41d480608fe0264b59f2861babaf7050a722` |
| Infernal | `container_infernal` | `quay.io/biocontainers/infernal@sha256:05ae1ca6cc76c27180524bc38c5b1e17adf9377be5b8c644d3e8e707848d4d99` |
| InterProScan | `container_interproscan` | `docker.io/interpro/interproscan@sha256:dc58b7c147fbbf00c2dd4f5ced42121fc1e8841fcbc7cc2c484380248ff76d11` |
| Liftoff | `container_liftoff` | `quay.io/biocontainers/liftoff@sha256:460d5e82b0c59e8348633f3e0b9a19cf29f9227f7457e90bd7f1d1a2403b3555` |
| Minimap2 + samtools | `container_minimap2_samtools` | `avelt/minimap2_samtools@sha256:70dcb87bb8021c90fc5eb660bbe1e6fc6bedadbf85c552c66704d27957b1f4ba` |
| Mikado | `container_mikado` | `quay.io/biocontainers/mikado@sha256:dd6f5a2a2d7fdbab73c835cd0f49bd1444ecaddf8e4cd96fbf0fe24f5ecf5f22` |
| MultiQC | `container_multiqc` | `quay.io/biocontainers/multiqc@sha256:b65e3fe879df27b92334dda0fd987a6e21bdee09a2848551d4f287099a93b7ac` |
| OMArk | `container_omark` | `quay.io/biocontainers/omark@sha256:84413cc19053c5d6452fbff245c9e6980b3f16aabdf991f9e51d7b9f2e0e0843` |
| PsiCLASS + samtools | `container_psiclass_samtools` | `avelt/psiclass_samtools@sha256:5cad8ecfd81293287bb6612ac8a6daaf17e626339016d326cf79615606acb285` |
| Python | `container_python` | `docker.io/library/python@sha256:57cd7c3a7a273101a6485ba99423ee568157882804b1124b4dd04266317710de` |
| Salmon | `container_salmon` | `quay.io/biocontainers/salmon@sha256:71ffc3b4961971159a6a2327d55686fb499c43335644ea5623476a082e826fc0` |
| SQANTI3 | `container_sqanti3` | `quay.io/biocontainers/sqanti3@sha256:3bd6ec96b3f1c9cae69cfef54ba0522b7d99efa7ebb0ff6a611841aa6784f74c` |
| STAR | `container_star` | `quay.io/biocontainers/star@sha256:f5910f39a9f5bc171a51fe7400d33e7586cb353c47d759a7c190562322150067` |
| StringTie | `container_stringtie` | `avelt/stringtie@sha256:856395c26e0c36544ef5c66e24badcac4f68fd5fa51864a0f964a737250545bb` |
| TransDecoder | `container_transdecoder` | `quay.io/biocontainers/transdecoder@sha256:c70f3a30cc8f3aecccb1d8978b9a49865d3994ebd7885361ab4c9dd820bd17f5` |
| tRNAscan-SE | `container_trnascan` | `quay.io/biocontainers/trnascan-se@sha256:e573090368974ff1228e6894828c6c8a132dfecc3198f5e9fb76832f8f434f29` |
| CPAT | `container_cpat` | `quay.io/biocontainers/cpat@sha256:87366fff67d441f64e0ac4681ccbaf1147f2c0601f3df86bb99f228d7f9a9000` |
| FLAIR | `container_flair` | `quay.io/biocontainers/flair@sha256:187e2e22535d73ecc724afc7e474d9908b6e43a55f8588e8566db3bea2eba79e` |
| BUSCO | `container_busco` | `quay.io/biocontainers/busco@sha256:d55ad622a5cafcd63c42fc309108688ab255bb9586ee756a5149e249d418c8bd` |

## Dockerfile base images

Only Dockerfiles used to build TITAN-maintained `avelt/*` images are tracked.
They are not used by the default workflow at runtime, but their base images are
locked too:

| Original base | Locked base |
| --- | --- |
| `ubuntu:24.10` | `ubuntu@sha256:f995e05e8adc3292853cc37e6edda72351f8002ce7469a29322d19e01529cb9f` |
| `ubuntu:20.04` | `ubuntu@sha256:c664f8f86ed5a386b0a340d981b8f81714e21a8b9c73f658c4bea56aa179d54a` |
| `rockylinux:9` | `rockylinux@sha256:d644d203142cd5b54ad2a83a203e1dee68af2229f8fe32f52a30c6e1d3c3a9e0` |

## Test scope

The minimal dataset under `test-data/minimal` is realistic for contracts: assemblies, Liftoff annotation, RNA-seq single/paired/long reads, protein samplesheet, EGAPx YAML and AEGIS evidence fixtures.

Validate container pinning and profile resolution with:

```bash
python3 scripts/validate_container_pins.py
python3 scripts/validate_profiles.py
python3 scripts/validate_minimal_test_data.py
nextflow config -profile test
nextflow run main.nf -profile test -stub-run -ansi-log false
```

The `-stub-run` exercises every TITAN process and the full workflow graph on the realistic fixtures. It does not run the scientific tools themselves. A real full run of EDTA, BRAKER3, EGAPx and AEGIS requires production-scale databases, licenses/resources and enough container storage; it should be executed on the target HPC/container runtime after the same pin validation.
