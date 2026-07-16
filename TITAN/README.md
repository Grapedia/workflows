# TITAN: The Intensive Transcript Annotation Pipeline


**Contributors:**

David Navarro, Antonio Santiago, José Tomás Matus [TomsBio Lab](https://tomsbiolab.com/)

Amandine Velt, Camille Rustenholz [INRAE](https://eng-svqv.colmar.hub.inrae.fr/)

Marco Moretto [Fondazione Edmund Mach](https://fmach.it/)

_____________________________________________________________________________________________

📖 **Full documentation**: [Read the Docs](https://grapedia.readthedocs.io/en/latest/workflows.html#titan-the-intensive-transcript-annotation-pipeline)

🛠️ **Installation and input preparation**: [docs/user/installation.md](docs/user/installation.md)

_____________________________________________________________________________________________

## **Inputs**

The following parameters are defined in `nextflow.config` and are required for the pipeline execution.

### **Workflow's general parameters**
- **`output_dir`**: Path to output directory, where the final files will be write.

TITAN has one public execution contract: evidence generation and Aegis are always launched together in the same Nextflow graph. The former `generate_evidence_data`, `aegis` and `all` workflow modes are no longer user-facing options.

### **Genome Assemblies**
- **`previous_assembly`**: Path to the **reference genome assembly** (FASTA).  
  _Example_: `data/assemblies/T2T_ref.fasta`
- **`new_assembly`**: Path to the **new genome assembly** to be annotated (FASTA).  
  _Example_: `data/assemblies/riesling.hap1.chromosomes.phased.fa`

### **Annotations & Data evidences**
- **`previous_annotations`**: Path to the **GFF3 file** containing annotations for the previous assembly.  
  _Example_: `data/annotations/PN40024_5.1_on_T2T_ref_with_names.gff3`
- **`RNAseq_samplesheet`**: Path to the **RNA-seq samplesheet** (CSV file) listing RNA-seq datasets to be used.  
  _Example_: `data/RNAseq_data/RNAseq_samplesheet.csv`
- **`RNAseq_data_dir`**: Directory containing FASTQ/FASTA files referenced by RNA-seq `sampleID` values.
  _Example_: `data/RNAseq_data`
- **`protein_samplesheet`**: Path to the **protein data samplesheet** (CSV file) listing protein datasets to be used.  
  _Example_: `data/protein_data/samplesheet.csv`
- **`egapx_paramfile`**: Path to the **NCBI/EGAPx YAML parameter file**. EGAPx is mandatory in TITAN.
  _Example_: `data/input_egapx.yaml`

### **Other Parameters**
- **EDTA**: EDTA is not a runtime option anymore. TITAN always runs EDTA because Aegis requires a hard-masked genome.
- **EGAPx**: EGAPx is not a runtime option anymore. TITAN always runs EGAPx from `egapx_paramfile` and publishes its results under `${output_dir}/egapx`.
- **Long reads**: Long-read processing is not controlled by a flag anymore. TITAN detects it from `RNAseq_samplesheet` rows where `library_layout` is `long`.
- **`edta_cpus`**: CPU allocation for the EDTA process.
  default: `"5"`
- **`egapx_cpus`**: CPU allocation for the EGAPx process.
  default: `"5"`
- **`egapx_version` / `egapx_revision` / `egapx_container`**: EGAPx runner and container version. Defaults are `0.5.2`, `v0.5.2` and a digest-pinned `ncbi/egapx` image.
- **`aegis_version` / `aegis_container`**: AEGIS CLI and container version. Defaults are `v0.3.25` and a digest-pinned `tomsbiolab/aegis` image.
- **`PSICLASS_vd_option`**: For PSICLASS process, the minimum average coverage depth of a transcript to be reported (FLOAT). This option is used to reduce the number of false monoexon genes.
  default: `"5.0"`
- **`PSICLASS_c_option`**: For PSICLASS process, only use the subexons with classifier score <= than the given number (FLOAT). This option is used to reduce the number of false monoexon genes.
  default: `"0.03"`
- **`STAR_memory_per_job`**: For STAR alignment process. If the depth of your RNAseq samples is high, TITAN may crash with an out of memory error, during the STAR alignment step. You can increase the memory here, it's in bytes, for example 60000000000 is about 55Gb per sample/job.
  default: `"60000000000"`
### **Input contract inventory**

Current static inventory from `main.nf`, `workflows/titan.nf`, `subworkflows/*.nf` and `modules/*.nf`:

| Input | Parameter or channel | Expected content | Used by |
| --- | --- | --- | --- |
| Target assembly | `--new_assembly` | FASTA file for the genome to annotate | Liftoff, AGAT, STAR, HISAT2, Minimap2, EDTA, BRAKER3, Aegis |
| Reference assembly | `--previous_assembly` | FASTA file matching the previous annotation | Liftoff |
| Previous annotation | `--previous_annotations` | GFF3 annotation on the reference assembly | Liftoff, then AGAT CDS extraction |
| RNA-seq samplesheet | `--RNAseq_samplesheet` | CSV with header `sampleID,SRA_or_FASTQ,library_layout`; `library_layout` values are `single`, `paired` or `long` | read preparation, fastp, Salmon, STAR, HISAT2, Minimap2 |
| RNA-seq data directory | `--RNAseq_data_dir` | Directory containing FASTQ/FASTA files matching RNA-seq sample IDs | read preparation, fastp, Minimap2 |
| Protein samplesheet | `--protein_samplesheet` | CSV with header `organism,filename`; filenames may be relative paths in the minimal fixtures, while historical modules still mount `data/protein_data` | BRAKER3 |
| EGAPx parameter file | `--egapx_paramfile` | YAML parameter file for EGAPx | Mandatory EGAPx process |
| Tool/resource options | `--edta_cpus`, `--egapx_cpus`, `--PSICLASS_vd_option`, `--PSICLASS_c_option`, `--STAR_memory_per_job` | Tool resource and tuning values | Tool commands |

`library_layout` is split into `single`, `paired` and `long` branches. `single` and `paired` feed the short-read path; the presence of at least one `long` row automatically enables Minimap2/StringTie long-read processing and the long-read BRAKER3/Aegis branch. EDTA and EGAPx are always part of TITAN.

### **Named Evidence Contract**

The internal evidence-generation subworkflow emits named evidence channels instead of a mixed key/file channel. The Aegis subworkflow consumes explicit inputs for the hard-masked genome, Liftoff, EGAPx GFF3, BRAKER/AUGUSTUS, GeneMark, STAR/StringTie, STAR/PsiCLASS and detected long-read evidence.

Aegis consumes those channels directly. The EDTA hard-masked genome flows as `EDTA.masked_genome -> generate_evidence_data.masked_genome -> aegis(masked_genome)`, so TITAN no longer depends on `assembly_masked.EDTA.fasta` being rediscovered from `output_dir`.

### **Parameter validation**

TITAN validates required parameters at workflow startup, before building input channels or launching heavy processes. The following parameters must be provided and point to existing files where applicable:

```text
--output_dir
--egapx_paramfile
--RNAseq_samplesheet
--protein_samplesheet
--new_assembly
--previous_assembly
--previous_annotations
```

Invalid input fails early with explicit messages such as:

```text
Missing required parameter(s): --RNAseq_samplesheet
Required input file(s) not found:
  --previous_annotations: test-data/minimal/valid/missing.gff3
--workflow is no longer supported. TITAN always runs evidence generation followed by Aegis in one graph.
```

### **Output contract inventory**

| Output family | Main files or patterns | Output location |
| --- | --- | --- |
| Liftoff annotation | `liftoff_previous_annotations.gff3`, `unmapped_features.txt` | `${output_dir}` |
| CDS FASTA from annotation | `${genome}.CDS.fasta.gz` | `${output_dir}/intermediate_files/liftoff/gff3_to_cds_fasta` |
| Salmon index and strand inference | `salmon_index/`, `${sample_ID}.strand_info.classified` | `${output_dir}/intermediate_files/salmon_index`, `${output_dir}/intermediate_files/salmon_strand` |
| Prepared and trimmed reads | `*.trimmed.fastq.gz` | `${output_dir}/intermediate_files/evidence_data/RNAseq_data/trimmed_data` |
| Genome indices | STAR `${genome}_index`, HISAT2 `${genome}.*.ht2`, Minimap2 `${genome}.mmi` | `${output_dir}/intermediate_files/evidence_data/*_databases` |
| Alignments | STAR/HISAT2/Minimap2 sorted BAM files | `${output_dir}/intermediate_files/evidence_data/RNAseq_alignments` |
| Transcript assemblies | per-sample `*_transcriptome.gtf`, `*_transcriptome.AltCommands.gtf`, `${sample_ID}_vote.gtf` | `${output_dir}/intermediate_files/evidence_data/transcriptomes` and `${output_dir}/intermediate_files/transcriptomes` |
| Merged transcriptomes | `merged_transcriptomes.*.gtf`, `stranded_merged_output.combined.gtf`, optional `unstranded_merged_output.combined.gtf` | `${output_dir}`, `${output_dir}/tmp` and intermediate transcriptome directories |
| EDTA | `*TElib.fa`, `*TEanno.gff3`, `*MAKER.masked` | `${output_dir}/tmp` |
| BRAKER3 | `augustus.hints.gff3`, `genemark.gtf`, `genemark_supported.gtf`, `braker.gff3` | `${output_dir}` |
| Aegis final annotation | `final_annotation.gff3`, `final_annotation_proteins_all.fasta`, `final_annotation_proteins_main.fasta` | `${output_dir}/aegis_outputs` |
| Diamond2GO | `*-diamond*` | `${output_dir}/Diamond2GO_outputs` |
| EGAPx | `egapx.complete.genomic.gff3`, `egapx.complete.genomic.gtf`, `egapx.complete.proteins.faa`, `egapx.complete.cds.fna`, `egapx.complete.transcripts.fna`, `egapx.annotated_genome.asn`, `egapx_out/`, `versions.yml` | `${output_dir}/egapx` |
| Provenance | `evidence_manifest.json`, `versions.yml` | `${output_dir}/provenance` |

## **Example command-line to run**

```bash
#!/usr/bin/env bash
# Navigate to the project workflow directory
cd /path/to/projectDir/workflows/TITAN
# Load required Nextflow module
module load nextflow/24.04.3
# Run the full workflow and generate its DAG
nextflow run main.nf \
  -with-dag dag_titan.png
```

## **Local test profile**

TITAN includes a minimal local test profile for configuration and lightweight workflow checks without Slurm, Docker or production data.

```bash
nextflow config -profile test
python3 scripts/validate_profiles.py
nextflow run main.nf -profile test -stub-run -ansi-log false
```

The `test` profile uses synthetic fixtures under `test-data/minimal/valid`, `RNAseq_data_dir = test-data/minimal/valid/rnaseq`, and ignored transient directories `test-results/` and `test-work/`. `-stub-run` exercises the full graph, including mandatory EDTA/EGAPx, auto-detected long reads and direct EDTA-to-Aegis channel wiring, without running scientific containers.

Production profiles are separated from the test profile:

```bash
nextflow config -profile local
nextflow config -profile apptainer
nextflow config -profile slurm,apptainer
```

Use `local` for a local Docker run, `apptainer` for a local Apptainer run, and `slurm,apptainer` for HPC production. For test-data resolution checks with HPC settings, use `test,slurm,apptainer` so the test input paths are applied before the Slurm/Apptainer runtime overrides.

Validate the synthetic fixture set directly with:

```bash
python3 scripts/validate_minimal_test_data.py
sha256sum -c test-data/minimal/checksums.sha256
```

Development notes for the completed P0 hardening work are in `docs/development/p0-hardening.md`; the detailed audit and inventory are in `docs/development/audit.md`; the architecture audit and refactor target are in `docs/development/architecture-audit.md`.

## Workflow DAG

![Workflow Diagram](data_example/TITAN_diagram.jpg)

## Workflow Components

### **Liftoff Annotations**
Generates a **GFF3 file**, which is used by **Aegis** in the final step.

### **NCBI/egapx Annotations**
Runs the **NCBI EGAPx** workflow from `egapx_paramfile` using the official EGAPx runner `v0.5.2` and the digest-pinned official `ncbi/egapx` Docker image. TITAN publishes named EGAPx outputs under `${output_dir}/egapx`: GFF3, GTF, proteins, CDS, transcripts, ASN, full EGAPx output directory and `versions.yml`. The EGAPx GFF3 is passed to AEGIS as an additional annotation evidence.

### **StringTie Merging (Short Reads - HISAT2)**
Generates **GTF file(s)**.  
- For each **short-read RNA-seq dataset**, a transcriptome is assembled using **HISAT2/StringTie**.  
- Transcriptomes are then **merged separately** for **stranded** and **unstranded** samples (unstranded is optional).  
- ⚠️ **This output is not used yet**—it is intended for **lncRNA detection**, which is not yet implemented.

### **StringTie Merging (Short Reads - STAR)**
Generates **GTF file(s)**, which are used by **Aegis** in the final step.  
- For each **short-read RNA-seq dataset**, a transcriptome is assembled using **STAR/StringTie** with **default** and **alternative (alt) parameters**.  
- Transcriptomes are **merged separately** for **stranded** and **unstranded** samples (unstranded is optional) and for **default and alt parameters** separately.

### **StringTie Merging (Long Reads)**
Generates a **single GTF file**, which is used by **Aegis** in the final step.  
- For each **IsoSeq RNA-seq dataset**, a transcriptome is assembled using **Minimap2/StringTie** with **default** and **alt parameters**.  
- These transcriptomes are then **merged into a single GTF file**, for **default and alt parameters** separately.
- Long-read integration is automatic when `RNAseq_samplesheet` contains at least one row where `library_layout` is `long`.

### **BRAKER3 Gene Prediction**
#### **braker3_prediction OR braker3_prediction_with_long_reads**
- **BRAKER3** is run with all **BAM files**, generated using **STAR (short reads)** or **Minimap2 (long reads)**.
- It also incorporates **protein FASTA files** to generate:
  - A **GeneMark GTF file**.
  - An **AUGUSTUS GFF3 file**, both of which are used by **Aegis** in the final step.

### **GFFCompare**
Generates a **GTF file**, which is used by **Aegis** in the final step.  
- For each **short-read RNA-seq dataset**, a transcriptome is assembled using **STAR/PsiCLASS**.  
- Transcriptomes are **merged separately** for **stranded** and **unstranded** samples using **GFFCompare** (unstranded is optional).

### **EDTA**
*(Integrated into TITAN, though not explicitly shown in the DAG)*  
Generates a **FASTA file** of the **hard-masked genome assembly**, which is used by **Aegis** in the final step.  

- **Input**: The **FASTA file** of the genome assembly to be annotated.  
- **Output**: A **hard-masked version** of the genome assembly.  

---

## **Final Integration by Aegis**
At the final step, **AEGIS** runs `aegis merge` on the named annotation evidence, then `aegis extract` to derive the final protein FASTA files. Merge priority follows the command-line input order: Liftoff, AUGUSTUS, GeneMark, EGAPx, STAR/StringTie, long-read StringTie when present, and STAR/PsiCLASS evidence.

AEGIS integrates multiple annotation sources:
✔ **Hard-masked genome assembly**
✔ **Liftoff annotations**  
✔ **NCBI/egapx annotations**  
✔ **GeneMark annotations**  
✔ **AUGUSTUS annotations**  
✔ **Stranded (and unstranded, if available) annotations from STAR/StringTie**  
✔ **IsoSeq annotations (if available)**  
✔ **Stranded (and unstranded, if available) annotations from STAR/PsiCLASS**  
### **📖 Reference**  
<https://github.com/Tomsbiolab/aegis>

This results in the final **GFF3 annotation file**.

## **Diamond2GO**  
*(Integrated into TITAN, though not explicitly shown in the DAG)*  

Diamond2GO performs **functional gene annotation** based on the **final Aegis protein FASTA files** (both "main" and "all" sets).  

### **🔹 Input**  
- FASTA files containing the **"main"** and **"all"** protein sequences generated by Aegis.  

### **🔹 Output**  
- Functional annotations assigned by **Diamond2GO**.  

## **📖 Tools version used**  

### Aegis
- Source from [GitHub](https://github.com/Tomsbiolab/aegis)
- **Version**: v0.3.25 from the current Docker image label
- **Docker image**: `tomsbiolab/aegis@sha256:de88470b3fb4fbab3ff2d5fa0fb9fed36b55952d1e383d3fdb2f5a3a530d84e6`

### Agat
- **Version**: 1.2.0
- **Docker image**: `quay.io/biocontainers/agat@sha256:7ea8fa5a8428758cd87e3a5dcfaf277febdfcae95cd1fe473770abf8b928ec99`

### BRAKER3
- **Version**: v3.0.8  
- **Dependencies**:
  - **spaln**: Ver.2.3.3  
  - **ProtHint**: 2.6.0  
  - **GeneMark-ETP**: Last version from [GitHub](https://github.com/gatech-genemark/GeneMark-ETP) (commit: `81ac83e`)  
  - **Augustus**: Last version from [GitHub](https://github.com/Gaius-Augustus/Augustus) (commit: `487b12b`)  
  - **bamtools**: v2.5.2  
  - **samtools**: 1.9  
  - **diamond**: v2.1.9  
  - **cdbfasta**: Last version from [GitHub](https://github.com/gpertea/cdbfasta.git) (commit: `da8f5ba`)  
  - **TSEBRA**: Last version from [GitHub](https://github.com/Gaius-Augustus/TSEBRA) (commit: `c87ba3a`)
- **Docker image**: `avelt/braker3@sha256:e69a9aaaafa81e4da5b2bbb98ae120d873018ae40453630f60051ecd5f622c44`

### Diamond2GO
- **Version**: Last version from [GitHub](https://github.com/rhysf/Diamond2GO.git) (commit: `57bb4cc`)
- **Docker image**: `avelt/diamond2go@sha256:40f1063307f98a2357d60b306bd7d79b6088591c1613e6552613da24002e8360`

### EDTA
- **Version**: Last version from [GitHub](https://github.com/oushujun/EDTA.git)  
- **Dependencies**:
  - **bedtools**: v2.30.0  
  - **samtools**: 1.9
- **Docker image**: `avelt/edta@sha256:607529be8e85c5b13dbed44135b35a0791ff9f8df41d4ea71169d47179315044`

### egapx
- **Version**: 0.5.2 [GitHub](https://github.com/ncbi/egapx)
- **Dependencies**:
  - **Python**: 3.11+ for the EGAPx runner
  - **Nextflow**: v23.10.1+ for the nested EGAPx workflow
- **Docker image**: `ncbi/egapx@sha256:bc657b232d93364d5f3b75ad3bfaf14b6267e46173672b609f26078d48a04298`

### fastp
- **Version**: 0.23.2
- **Docker image**: `quay.io/biocontainers/fastp@sha256:0bdf8d8254fc86dd9038551d68dbcb72562e65560b9ce0ea08c1329d2f8587b4`

### GFFCompare
- **Version**: 0.12.6
- **Docker image**: `avelt/gffcompare@sha256:bd411c13352a2545641c8c34b701030b3977056696607e85b4e86c876d10a82c`

### HISAT2
- **Version**: 2.2.1
- **Docker image**: `avelt/hisat2@sha256:022933fd0d30fe9fdfd83c175f7e41d480608fe0264b59f2861babaf7050a722`

### Liftoff
- **Version**: 1.5.1  
- **Docker image**: `quay.io/biocontainers/liftoff@sha256:460d5e82b0c59e8348633f3e0b9a19cf29f9227f7457e90bd7f1d1a2403b3555`

### Minimap2
- **Version**: 2.28  
- **Dependencies**:
  - **samtools**: 1.9
- **Docker image**: `avelt/minimap2_samtools@sha256:70dcb87bb8021c90fc5eb660bbe1e6fc6bedadbf85c552c66704d27957b1f4ba`

### PsiCLASS
- **Version**: 1.0.2  
- **Docker image**: `avelt/psiclass_samtools@sha256:5cad8ecfd81293287bb6612ac8a6daaf17e626339016d326cf79615606acb285`

### sra-tools
- **Version**: 3.1.1  
- **Docker image**: `quay.io/biocontainers/sra-tools@sha256:05de2c580cccc4c609ec7c645902563e5d5ffbd366662e1983cb152545ec7bc0`

### STAR
- **Version**: 2.7.11b
- **Docker image**: `quay.io/biocontainers/star@sha256:f5910f39a9f5bc171a51fe7400d33e7586cb353c47d759a7c190562322150067`

### Salmon
- **Version**: 1.10.3
- **Docker image**: `quay.io/biocontainers/salmon@sha256:71ffc3b4961971159a6a2327d55686fb499c43335644ea5623476a082e826fc0`

### StringTie
- **Version**: 2.2.3  
- **Docker image**: `avelt/stringtie@sha256:856395c26e0c36544ef5c66e24badcac4f68fd5fa51864a0f964a737250545bb`

## **📖 Reference**  

**Aegis** : TO DO

**Agat** : Dainat J. *Another Gtf/Gff Analysis Toolkit (AGAT): Resolve interoperability issues and accomplish more with your annotations.*  
📄 [Plant and Animal Genome XXIX Conference, 2022](https://github.com/NBISweden/AGAT)

**BRAKER3** : Gabriel et al. *BRAKER3: Fully automated genome annotation using RNA-seq and protein evidence with GeneMark-ETP, AUGUSTUS, and TSEBRA*  
📄 [Genome Res, 2024](https://pubmed.ncbi.nlm.nih.gov/38866550/)

**Diamond2GO** : Golden et al. *DIAMOND2GO: A rapid Gene Ontology assignment and enrichment tool for functional genomics.*  
📄 [bioRxiv, 2024](https://www.biorxiv.org/content/10.1101/2024.08.19.608700v1)

**EDTA** : Ou et al. *Benchmarking Transposable Element Annotation Methods for Creation of a Streamlined, Comprehensive Pipeline.*  
📄 [Genome Biol., 2019](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1905-y)

**egapx (NCBI)** : *The NCBI Eukaryotic Genome Annotation Pipeline.*  
📄 [ncbi.nlm.nih.gov](https://www.ncbi.nlm.nih.gov/refseq/annotation_euk/process)

**fastp** : Chen et al. *fastp: an ultra-fast all-in-one FASTQ preprocessor*  
📄 [Bioinformatics., 2018](https://academic.oup.com/bioinformatics/article/34/17/i884/5093234)

**GFFCompare** : Pertea et al. *GFF Utilities: GffRead and GffCompare.*  
📄 [F1000Research., 2020](https://pubmed.ncbi.nlm.nih.gov/32489650/)

**GffRead** : Pertea et al. *GFF Utilities: GffRead and GffCompare.*  
📄 [F1000Research., 2020](https://pubmed.ncbi.nlm.nih.gov/32489650/)

**HISAT2** : Kim et al. *Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype*  
📄 [Nature biotechnology, 2019](https://www.nature.com/articles/s41587-019-0201-4)

**Liftoff** : Shumate et al. *Liftoff: accurate mapping of gene annotations*  
📄 [Bioinformatics, 2021](https://academic.oup.com/bioinformatics/article/37/12/1639/6035128?login=false)

**Minimap2** : Li et al. *Minimap2: pairwise alignment for nucleotide sequences*  
📄 [Bioinformatics, 2018](https://academic.oup.com/bioinformatics/article/34/18/3094/4994778)

**PsiCLASS** : Song et al. *A multi-sample approach increases the accuracy of transcript assembly*  
📄 [Nature communications., 2019](https://www.nature.com/articles/s41467-019-12990-0)

**sra-tools** : [GitHub repository](https://github.com/ncbi/sra-tools)

**STAR** : Dobin et al. *STAR: ultrafast universal RNA-seq aligner*  
📄 [Bioinformatics, 2013](https://academic.oup.com/bioinformatics/article/29/1/15/272537)

**Salmon** : Patro et al. *Salmon provides fast and bias-aware quantification of transcript expression*  
📄 [Nature methods, 2017](https://www.nature.com/articles/nmeth.4197)

**Samtools** : Li et al. *The Sequence Alignment/Map format and SAMtools*  
📄 [Bioinformatics, 2009](https://pmc.ncbi.nlm.nih.gov/articles/PMC2723002/)

**StringTie** : Shumate et al. *Improved transcriptome assembly using a hybrid of long and short reads with StringTie*  
📄 [PLOS Computational Biology, 2022](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009730)
