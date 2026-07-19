# TITAN Tool Commands Reference

## Table Of Contents

- [Execution Order](#execution-order)
- [Input Validation](#input-validation)
- [tRNA Annotation](#trna-annotation)
- [Rfam Annotation](#rfam-annotation)
- [Short-Read Staging](#short-read-staging)
- [Liftoff Transfer](#liftoff-transfer)
- [AGAT Conversion For Liftoff Evidence](#agat-conversion-for-liftoff-evidence)
- [EGAPx](#egapx)
- [Salmon Strand Inference](#salmon-strand-inference)
- [STAR Index And Alignment](#star-index-and-alignment)
- [HISAT2 Index And Alignment](#hisat2-index-and-alignment)
- [Long-Read Staging And Minimap2](#long-read-staging-and-minimap2)
- [StringTie Transcript Assembly](#stringtie-transcript-assembly)
- [PsiCLASS Transcript Assembly](#psiclass-transcript-assembly)
- [Merging Transcript Assemblies](#merging-transcript-assemblies)
- [GffCompare For PsiCLASS](#gffcompare-for-psiclass)
- [FLAIR Isoforms](#flair-isoforms)
- [EDTA Repeat Annotation](#edta-repeat-annotation)
- [Protein FASTA Normalization](#protein-fasta-normalization)
- [BRAKER3](#braker3)
- [Helixer](#helixer)
- [AEGIS Evidence Merge](#aegis-evidence-merge)
- [Diamond2GO](#diamond2go)
- [eggNOG-mapper](#eggnog-mapper)
- [InterProScan](#interproscan)
- [Mikado And TransDecoder](#mikado-and-transdecoder)
- [lncRNA Candidate Annotation](#lncrna-candidate-annotation)
- [SQANTI3 Long-Read QC](#sqanti3-long-read-qc)
- [Final Annotation Validation](#final-annotation-validation)
- [Final Transcriptome And Expression Support](#final-transcriptome-and-expression-support)
- [BUSCO](#busco)
- [OMArk](#omark)
- [AGAT Structural Statistics](#agat-structural-statistics)
- [ncRNA QC](#ncrna-qc)
- [MultiQC](#multiqc)
- [Provenance Manifests](#provenance-manifests)
- [Skipped And Sentinel Processes](#skipped-and-sentinel-processes)
- [Process Coverage Inventory](#process-coverage-inventory)


This reference documents the commands launched by TITAN in logical workflow
order. It is a technical companion to `docs/reference/tools.md`: that file
explains what each tool contributes, while this file records the command lines,
important options and branch behavior.

Command snippets use placeholders such as `<genome>`, `<sample_ID>`,
`<cpus>` and `<params.*>` for Nextflow inputs and runtime parameters. Optional
branches emit empty sentinel outputs when disabled so downstream aggregation
remains stable.

## Execution Order

The order below follows `workflows/titan.nf`,
`subworkflows/generate_evidence_data.nf` and `subworkflows/aegis.nf`.

1. Validate inputs and stage optional reference data.
2. Run early ncRNA branches on the target genome: tRNAscan-SE and Rfam.
3. Prepare RNA-seq reads, transfer the previous annotation and run EGAPx.
4. Build RNA-seq indices, infer strandedness, align reads and assemble
   transcript evidence.
5. Run EDTA, normalize protein evidence and run BRAKER3.
6. Run optional Helixer.
7. Merge evidence with AEGIS and extract final proteins.
8. Run functional annotation on AEGIS proteins.
9. Run optional Mikado/TransDecoder, lncRNA and SQANTI3 branches.
10. Validate and summarize the final annotation with expression, BUSCO, OMArk,
    AGAT, ncRNA QC, MultiQC and provenance manifests.

## Input Validation

Process: `validate_inputs`  
Module: `modules/validate_inputs.nf`

When `--run_lncrna true` and `--lncrna_require_cpat_model true`, TITAN first
ensures the CPAT plant model exists, downloading it if needed:

```bash
bash scripts/download_cpat_plant_lncpipe.sh \
  --model-dir <params.cpat_model_dir>
```

Then it validates the complete user-facing input contract:

```bash
python3 scripts/validate_inputs.py \
  --project-dir <projectDir> \
  --new-assembly <new_assembly> \
  --previous-assembly <previous_assembly> \
  --previous-annotations <previous_annotations> \
  --rnaseq-samplesheet <RNAseq_samplesheet> \
  --rnaseq-data-dir <RNAseq_data_dir> \
  --protein-samplesheet <protein_samplesheet> \
  --egapx-paramfile <egapx_paramfile> \
  --egapx-executor <params.egapx_executor> \
  --psiclass-vd <params.PSICLASS_vd_option> \
  --psiclass-c <params.PSICLASS_c_option> \
  --run-eggnog-mapper <params.run_eggnog_mapper> \
  --eggnog-data-dir <params.eggnog_data_dir> \
  --run-helixer <params.run_helixer> \
  --helixer-model-dir <params.helixer_model_dir> \
  --helixer-model <params.helixer_model> \
  --run-interproscan <params.run_interproscan> \
  --interproscan-data-dir <params.interproscan_data_dir>
```

Important options:

| Option | Purpose | Why it matters |
|---|---|---|
| `--new-assembly`, `--previous-assembly`, `--previous-annotations` | Validate the core comparative annotation inputs. | Prevents late failures in Liftoff, EDTA and AEGIS. |
| `--rnaseq-samplesheet`, `--rnaseq-data-dir` | Validate RNA-seq metadata and local read paths. | Determines short/long branches and read staging. |
| `--protein-samplesheet` | Validate external protein evidence. | BRAKER3 requires at least one protein FASTA. |
| `--egapx-paramfile`, `--egapx-executor` | Validate the nested EGAPx run configuration. | EGAPx runs outside normal containers and fails late without these checks. |
| Optional tool flags and data dirs | Validate enabled optional branches. | Keeps disabled branches lightweight while failing fast for enabled ones. |

## tRNA Annotation

Process: `trnascan_se`  
Module: `modules/trnascan_se.nf`

Runs when `--run_trnascan true`; otherwise TITAN writes empty tRNA outputs.

```bash
export TMPDIR=/tmp

tRNAscan-SE -E \
  -o trnascan.out \
  -f trnascan.struct \
  -s trnascan.isotype \
  -m trnascan.stats \
  --thread <cpus> \
  <genome>
```

Important options:

| Option | Purpose | Why it matters |
|---|---|---|
| `-E` | Eukaryotic search mode. | Matches TITAN's eukaryotic genome scope. |
| `-o`, `-f`, `-s`, `-m` | Write raw table, structure, isotype and stats outputs. | Feeds GFF3 conversion and ncRNA QC. |
| `--thread` | Parallel execution. | Uses the Nextflow task allocation. |

Process: `trnascan_to_gff3`

```bash
python3 scripts/trnascan_to_gff3.py trnascan.out > trna.gff3
```

This internal converter turns the tRNAscan-SE table into a compact GFF3 used
by lncRNA filtering and ncRNA QC.

## Rfam Annotation

Process: `rfam_split_genome`  
Module: `modules/infernal_rfam.nf`

```bash
python3 - <<'PY'
# Split <genome> into one FASTA file per sequence under rfam_genome_parts/.
PY
```

The split allows one Infernal task per chromosome/contig.

Rfam reference data source used for TITAN validation:

* Dataset: global Rfam covariance-model library, not a `Vitis vinifera`
  species-specific resource.
* Files required by TITAN: `Rfam.cm`, `Rfam.clanin` and the `cmpress` index
  files generated from `Rfam.cm`.
* Source URL used by the TITAN downloader and Vitis vinifera validation runs:
  <https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT>
* Direct `Rfam.clanin` URL:
  <https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.clanin>
* Preparation script:
  [`scripts/download_rfam_data.sh`](../../scripts/download_rfam_data.sh)

The same staged directory is then passed to TITAN with `--rfam_data_dir`; the
Vitis vinifera run therefore uses `Rfam.clanin` from the EBI Rfam `CURRENT`
release selected at download time.

```bash
scripts/download_rfam_data.sh \
  --data-dir <rfam_data_dir> \
  --container <params.container_infernal>
```

Process: `infernal_rfam_search`

Runs when `--run_rfam true`; otherwise TITAN writes empty table/log files.

```bash
cmsearch --cpu <cpus> \
  --tblout <sequence_id>.rfam_hits.tbl \
  --cut_ga \
  --rfam \
  --nohmmonly \
  --clanin <params.rfam_data_dir>/Rfam.clanin \
  <params.rfam_data_dir>/Rfam.cm \
  <genome_part> > <sequence_id>.rfam_search.out
```

Important options:

| Option | Purpose | Why it matters |
|---|---|---|
| `--cpu` | Parallel Infernal search. | Controls per-sequence task speed. |
| `--tblout` | Machine-readable hit table. | Parsed into TITAN's Rfam GFF3. |
| `--cut_ga` | Use curated gathering thresholds. | Avoids arbitrary E-value cutoffs. |
| `--rfam` | Rfam-specific search behavior. | Keeps Infernal aligned with Rfam conventions. |
| `--nohmmonly` | Do not use HMM-only filters. | Keeps covariance-model sensitivity. |
| `--clanin` | Clan annotation file. | Handles clan-level overlap logic. |

Process: `infernal_rfam_merge`

```bash
cat rfam_search_results/*.rfam_hits.tbl > rfam_hits.tbl
cat rfam_search_results/*.rfam_search.out > rfam_search.out
python3 scripts/rfam_tblout_to_gff3.py rfam_hits.tbl > rfam_ncrna.gff3
```

The merge keeps one final Rfam table, one combined log and one GFF3 track.

## Short-Read Staging

Process: `prepare_RNAseq_fastq_files_short`  
Module: `modules/prepare_RNAseq_fastq_files_short.nf`

For ENA/SRA paired-end reads:

```bash
python3 scripts/download_sra_fastq.py \
  --accession <sample_ID> \
  --layout paired \
  --outdir . \
  --timeout-seconds <params.ena_download_timeout_seconds> \
  --max-attempts <params.ena_max_download_attempts> \
  --retry-wait-seconds <params.ena_retry_wait_seconds> \
  --verify-md5
```

For ENA/SRA single-end reads, `--layout single` is used. If
`--ena_verify_md5 false`, TITAN uses `--no-verify-md5`. For local FASTQ rows,
the module checks expected files and copies them to `prepared_1.fastq.gz` and
`prepared_2.fastq.gz`; single-end rows get an empty gzipped second mate.

Important options:

| Option | Purpose | Why it matters |
|---|---|---|
| `--accession` | ENA/SRA run accession. | Resolves remote FASTQ files. |
| `--layout` | `paired`, `single` or `long`. | Controls expected file count and naming. |
| `--timeout-seconds`, `--max-attempts`, `--retry-wait-seconds` | Download resilience. | Reduces transient ENA/network failures. |
| `--verify-md5` / `--no-verify-md5` | Integrity checking mode. | Enables strict downloads when checksums are available. |

Process: `trimming_fastq`  
Module: `modules/trimming_fastq.nf`

Paired-end:

```bash
fastp --thread <cpus> \
  -i <read_1> \
  -I <read_2> \
  -o <sample_ID>_1.trimmed.fastq.gz \
  -O <sample_ID>_2.trimmed.fastq.gz \
  --json <sample_ID>.fastp.json \
  --html <sample_ID>.fastp.html
```

Single-end:

```bash
fastp --thread <cpus> \
  -i <read_1> \
  -o <sample_ID>_1.trimmed.fastq.gz \
  --json <sample_ID>.fastp.json \
  --html <sample_ID>.fastp.html
```

Single-end rows also receive an empty gzipped `<sample_ID>_2.trimmed.fastq.gz`
placeholder for downstream tuple stability.

## Liftoff Transfer

Process: `liftoff_annotations`  
Module: `modules/liftoff_annotations.nf`

```bash
liftoff \
  -g <previous_annotations_gff3> \
  -o liftoff_previous_annotations.gff3 \
  -u unmapped_features.txt \
  <new_assembly_fasta> \
  <previous_assembly_fasta>
```

Important options:

| Option | Purpose | Why it matters |
|---|---|---|
| `-g` | Previous GFF3 annotation. | Defines features to transfer. |
| `-o` | Transferred GFF3 output. | Becomes evidence for AEGIS and Mikado. |
| `-u` | Unmapped feature report. | Preserves transfer failures for audit. |

## AGAT Conversion For Liftoff Evidence

Process: `agat_convert_gff3_to_gtf`  
Module: `modules/agat_convert_gff3_to_gtf.nf`

```bash
awk -F'\t' 'NF == 9 || /^#/' <liftoff_gff3> > sanitized.gff3
agat_convert_sp_gff2gtf.pl \
  --gff sanitized.gff3 \
  -o <liftoff_basename>.gtf
```

The `awk` filter keeps AGAT from failing on malformed non-GFF records. The GTF
is used as splice-junction correction evidence for FLAIR.

Process: `clean_liftoff_gff3_for_agat`  
Module: `modules/clean_liftoff_gff3_for_agat.nf`

```bash
python3 scripts/clean_liftoff_gff3_for_agat.py \
  --input <liftoff_gff3> \
  --output cleaned.OK.gff3 \
  --removed-ids removed_feature_ids.txt
```

This internal cleaner normalizes Liftoff records before AGAT sequence
extraction.

Process: `agat_convert_gff3_to_cds_fasta`  
Module: `modules/agat_convert_gff3_to_cds_fasta.nf`

```bash
fold -w 80 <genome> > reformatted.fa
agat_sp_extract_sequences.pl \
  -g <cleaned_liftoff_gff3> \
  -f reformatted.fa \
  -o <genome>.CDS.fasta
gzip <genome>.CDS.fasta
```

The CDS FASTA seeds Salmon indexing for strand inference.

## EGAPx

Process: `egapx`  
Module: `modules/egapx.nf`

EGAPx is mandatory. TITAN launches the official nested EGAPx runner outside a
normal process container, because it must control the host Docker/Apptainer or
Nextflow executor.

If `--egapx_runner_dir` is provided:

```bash
cp -R <params.egapx_runner_dir>/. egapx_runner/
```

Otherwise TITAN downloads the requested runner revision:

```bash
curl -fsSL https://github.com/ncbi/egapx/archive/refs/tags/<params.egapx_revision>.tar.gz \
  | tar -xz --strip-components=1 -C egapx_runner
```

Then it pins the nested EGAPx process container and runs:

```bash
printf "process.container = '%s'\n" "<params.container_egapx>" \
  > egapx_runner/ui/assets/config/docker_image.config

<params.egapx_python> egapx_runner/ui/egapx.py \
  <egapx_paramfile> \
  -e <params.egapx_executor> \
  -w <workdir>/egapx_work \
  -o <workdir>/egapx_out \
  -lc <egapx_local_cache> \
  -dv <params.egapx_data_version> \
  -c <params.egapx_config_dir>
```

The `-c` option is added only when `--egapx_config_dir` is set.

Important options:

| Option | Purpose | Why it matters |
|---|---|---|
| `-e` | Nested executor command. | Lets EGAPx run under the configured host executor. |
| `-w` | EGAPx work directory. | Isolates nested Nextflow state inside the TITAN task. |
| `-o` | EGAPx output directory. | Standard location for official EGAPx products. |
| `-lc` | Local cache. | Avoids repeated downloads and supports offline staging. |
| `-dv` | EGAPx data version. | Reproducibility of EGAPx reference data. |
| `-c` | Extra EGAPx config directory. | Site-specific executor or runtime config. |

TITAN requires and copies `complete.genomic.gff`, `complete.genomic.gtf`,
`complete.proteins.faa`, `complete.cds.fna`, `complete.transcripts.fna` and
`annotated_genome.asn`.

## Salmon Strand Inference

Process: `salmon_index`  
Module: `modules/salmon_index.nf`

```bash
salmon index \
  -t <liftoff_cds_fasta.gz> \
  -i salmon_index
```

Process: `salmon_strand_inference`  
Module: `modules/salmon_strand_inference.nf`

Paired-end:

```bash
salmon quant \
  -i <salmon_index> \
  -l A \
  -p <cpus> \
  -1 <read_1> \
  -2 <read_2> \
  -o <sample_ID>_quant \
  --validateMappings \
  --skipQuant \
  2> <sample_ID>.log
```

Single-end:

```bash
salmon quant \
  -i <salmon_index> \
  -l A \
  -p <cpus> \
  -r <read_1> \
  -o <sample_ID>_quant \
  --validateMappings \
  --skipQuant \
  2> <sample_ID>.log
```

Important options:

| Option | Purpose | Why it matters |
|---|---|---|
| `-l A` | Autodetect library type. | Infers strandness from each sample. |
| `--validateMappings` | Selective alignment validation. | Improves mapping reliability for inference. |
| `--skipQuant` | Do not produce full quantification. | Uses Salmon only to infer library type quickly. |

TITAN maps Salmon library types to downstream classes:

| Salmon type | TITAN class |
|---|---|
| `IU`, `U` | `unstranded` |
| `ISF`, `SF`, `FR` | `stranded_forward` |
| `ISR`, `SR`, `RF` | `stranded_reverse` |
| anything else | `unstranded` fallback |

## STAR Index And Alignment

Process: `star_genome_indices`  
Module: `modules/star_genome_indices.nf`

```bash
STAR \
  --runThreadN <cpus> \
  --runMode genomeGenerate \
  --genomeDir <genome>_index \
  --genomeFastaFiles <genome_fasta> \
  --genomeSAindexNbases <STAR_genomeSAindexNbases> \
  --sjdbGTFfile <params.STAR_sjdbGTFfile>
```

`--sjdbGTFfile` is added only when `--STAR_sjdbGTFfile` is provided. If
`--STAR_genomeSAindexNbases false`, TITAN derives it from genome length and
caps it between 1 and 14.

Process: `star_alignment`  
Module: `modules/star_alignment.nf`

Paired-end:

```bash
STAR \
  --readFilesCommand zcat \
  --genomeDir <star_database> \
  --runThreadN <cpus> \
  --readFilesIn <read_1> <read_2> \
  --outFileNamePrefix <sample_ID>_ \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMunmapped Within \
  --outSAMattributes Standard \
  --outSAMstrandField intronMotif \
  --limitBAMsortRAM <task.memory.bytes>
```

Single-end uses the same command with only `<read_1>`.

Important options:

| Option | Purpose | Why it matters |
|---|---|---|
| `--readFilesCommand zcat` | Read gzipped FASTQ. | Avoids manual decompression. |
| `--outSAMtype BAM SortedByCoordinate` | Direct sorted BAM output. | Feeds StringTie, PsiCLASS and BRAKER3. |
| `--outSAMstrandField intronMotif` | Add strand field from splice motifs. | Helps downstream transcript assembly. |
| `--limitBAMsortRAM` | STAR sort memory limit. | Uses Nextflow memory allocation explicitly. |

## HISAT2 Index And Alignment

Process: `hisat2_genome_indices`  
Module: `modules/hisat2_genome_indices.nf`

```bash
/hisat2-2.2.1/hisat2-build \
  -p <cpus> \
  <genome_fasta> \
  hisat2_index/<genome>
```

Process: `hisat2_alignment`  
Module: `modules/hisat2_alignment.nf`

Paired-end:

```bash
/hisat2-2.2.1/hisat2 \
  -p <cpus> \
  -x <hisat2_index_prefix> \
  <strandness_option> \
  -1 <read_1> \
  -2 <read_2> \
  | samtools sort -o <sample_ID>_Aligned.sort.bam -

samtools index -@ <cpus> <sample_ID>_Aligned.sort.bam
```

Single-end uses `-U <read_1>`.

Strandness options:

| TITAN class | Paired-end | Single-end |
|---|---|---|
| `unstranded` | no option | no option |
| `stranded_forward` | `--rna-strandness FR` | `--rna-strandness F` |
| `stranded_reverse` | `--rna-strandness RF` | `--rna-strandness R` |

## Long-Read Staging And Minimap2

Process: `prepare_RNAseq_fastq_files_long`  
Module: `modules/prepare_RNAseq_fastq_files_long.nf`

For ENA/SRA long reads:

```bash
python3 scripts/download_sra_fastq.py \
  --accession <sample_ID> \
  --layout long \
  --outdir . \
  --timeout-seconds <params.ena_download_timeout_seconds> \
  --max-attempts <params.ena_max_download_attempts> \
  --retry-wait-seconds <params.ena_retry_wait_seconds> \
  --verify-md5
```

Local long-read rows can be `FASTQ` (`<sample_ID>.fastq.gz`) or `FASTA`
(`<sample_ID>.fasta`). TITAN records `READ_FORMAT=fastq` or `READ_FORMAT=fasta`
for Minimap2.

Process: `minimap2_genome_indices`  
Module: `modules/minimap2_genome_indices.nf`

```bash
minimap2 -d <genome>.mmi <genome_fasta>
```

Process: `minimap2_alignment`  
Module: `modules/minimap2_alignment.nf`

```bash
minimap2 \
  -t <cpus> \
  -ax splice:hq \
  -uf \
  <genome>.mmi \
  <long_reads.fastq_or_fasta> \
  | samtools sort -o <sample_ID>_Aligned.sorted.bam -
```

Important options:

| Option | Purpose | Why it matters |
|---|---|---|
| `-d` | Build a Minimap2 index. | Reused across long-read samples. |
| `-ax splice:hq` | Splice-aware high-quality transcript preset. | Appropriate for long-read RNA alignments. |
| `-uf` | Transcript-strand/long-read mapping mode. | Keeps spliced long-read alignment behavior. |
| `samtools sort` | Coordinate-sort BAM output. | Required by StringTie and BRAKER3. |

## StringTie Transcript Assembly

Processes:

* `assembly_transcriptome_star_stringtie`
* `assembly_transcriptome_hisat2_stringtie`
* `assembly_transcriptome_minimap2_stringtie`

Modules:

* `modules/assembly_transcriptome_star_stringtie.nf`
* `modules/assembly_transcriptome_hisat2_stringtie.nf`
* `modules/assembly_transcriptome_minimap2_stringtie.nf`

All three call the shared wrapper:

```bash
bash scripts/run_stringtie_transcriptome.sh \
  scripts/Stringtie.sh \
  scripts/Stringtie_AltCommands.sh \
  <cpus> \
  <bam_file> \
  <short_or_long> \
  <sample_ID>_transcriptome.gtf \
  <sample_ID>_transcriptome.AltCommands.gtf
```

Default short-read StringTie:

```bash
stringtie \
  -p <threads> \
  -o <output.gtf> \
  <bam>
```

Default long-read StringTie:

```bash
stringtie \
  -L \
  -p <threads> \
  -o <output.gtf> \
  <bam>
```

Alternative short-read StringTie:

```bash
stringtie \
  -f 0.99 \
  -m 120 \
  -a 15 \
  -j 3 \
  -c 3 \
  -s 4.75 \
  -g 50 \
  -t \
  -p <threads> \
  -o <output.gtf> \
  <bam>
```

Alternative long-read StringTie adds `-L`:

```bash
stringtie \
  -f 0.99 \
  -m 120 \
  -a 15 \
  -j 3 \
  -c 3 \
  -s 4.75 \
  -g 50 \
  -t \
  -L \
  -p <threads> \
  -o <output.gtf> \
  <bam>
```

Important options:

| Option | Purpose | Why it matters |
|---|---|---|
| `-L` | Long-read mode. | Uses StringTie behavior intended for long-read RNA alignments. |
| `-p` | Threads. | Matches the task CPU allocation. |
| `-f 0.99` | High isoform fraction threshold. | Produces a stricter alternative transcript set. |
| `-m 120` | Minimum transcript length. | Drops very short transcript fragments. |
| `-a 15` | Minimum junction anchor length. | Enforces stronger splice support. |
| `-j 3` | Minimum junction read coverage. | Reduces weak splice junctions. |
| `-c 3` | Minimum transcript coverage. | Reduces low-coverage transcripts. |
| `-s 4.75` | Single-exon minimum coverage. | Makes single-exon calls more conservative. |
| `-g 50` | Bundle gap. | Controls transcript clustering distance. |
| `-t` | Disable trimming at transcript ends. | Keeps the alternative profile reproducible. |

## PsiCLASS Transcript Assembly

Process: `assembly_transcriptome_star_psiclass`  
Module: `modules/assembly_transcriptome_star_psiclass.nf`

```bash
/PsiCLASS-1.0.2/psiclass \
  -p <cpus> \
  -b <star_sorted_bam> \
  -o <sample_ID> \
  --vd <params.PSICLASS_vd_option> \
  -c <params.PSICLASS_c_option> \
  --primaryParalog
```

If PsiCLASS exits non-zero and no trusted splice-junction files were produced,
TITAN emits an empty `<sample_ID>_vote.gtf`; otherwise the failure is fatal.

## Merging Transcript Assemblies

Processes:

* `Stringtie_merging_short_reads_STAR`
* `Stringtie_merging_short_reads_hisat2`
* `Stringtie_merging_long_reads`

Modules:

* `modules/Stringtie_merging_short_reads_STAR.nf`
* `modules/Stringtie_merging_short_reads_hisat2.nf`
* `modules/Stringtie_merging_long_reads.nf`

STAR and HISAT2 short-read merges are split into stranded and unstranded sets:

```bash
stringtie --merge \
  -o merged_transcriptomes.STAR.short_reads.default_args.stranded.gtf \
  stranded_default_gtfs.txt

stringtie --merge \
  -o merged_transcriptomes.STAR.short_reads.alt_args.stranded.gtf \
  stranded_alt_gtfs.txt
```

HISAT2 uses the same command with `merged_transcriptomes.hisat2...` output
names. Unstranded merges run only when non-empty unstranded GTFs exist.

Long-read merges:

```bash
stringtie --merge \
  -o merged_transcriptomes.minimap2.long_reads.default_args.gtf \
  default_gtfs.txt

stringtie --merge \
  -o merged_transcriptomes.minimap2.long_reads.alt_args.gtf \
  alt_gtfs.txt
```

If no long reads are present, `empty_long_read_evidence` copies empty sentinel
GTFs instead of running StringTie.

When long reads are present, these Minimap2/StringTie merged GTFs are the
primary long-read transcript evidence passed to AEGIS and Mikado. They are
produced regardless of `--run_flair`.

## GffCompare For PsiCLASS

Process: `gffcompare`  
Module: `modules/gffcompare.nf`

```bash
/gffcompare-0.12.6/gffcompare \
  -o stranded_merged_output \
  -i stranded_gtfs.txt
```

If unstranded PsiCLASS GTFs exist:

```bash
/gffcompare-0.12.6/gffcompare \
  -o unstranded_merged_output \
  -i unstranded_gtfs.txt
```

GffCompare produces merged PsiCLASS GTFs for AEGIS and Mikado.

## FLAIR Isoforms

Process: `flair_isoforms`  
Module: `modules/flair.nf`

Runs only when long reads are present and `--run_flair true`; otherwise
`flair_empty_evidence` emits empty sentinel GTF/FASTA files.

```bash
flair align \
  -g <genome> \
  -r <long_reads.fastq_or_fasta> \
  -o <sample_ID>.flair \
  --threads <cpus>

flair correct \
  -q <sample_ID>.flair.bed \
  -f <liftoff_gtf> \
  -o <sample_ID>.flair.corrected \
  --threads <cpus>

flair collapse \
  -g <genome> \
  -r <long_reads.fastq_or_fasta> \
  -q <sample_ID>.flair.corrected_all_corrected.bed \
  -o <sample_ID>.flair \
  --threads <cpus>
```

Process: `flair_merge_isoforms`

```bash
cat flair_gtfs/*.gtf > flair_isoforms.merged.gtf
cat flair_fastas/*.fa flair_fastas/*.fasta > flair_isoforms.merged.fasta
```

Important options:

| Option | Purpose | Why it matters |
|---|---|---|
| `flair align` | Align long-read isoforms. | Generates BED evidence for correction. |
| `flair correct -f <liftoff_gtf>` | Correct splice junctions against transferred annotation. | Uses prior annotation to improve isoform confidence. |
| `flair collapse` | Collapse corrected reads into isoforms. | Produces an additional optional GTF/FASTA evidence track passed to AEGIS and Mikado when `--run_flair true`. |

FLAIR does not replace the Minimap2/StringTie long-read evidence. If
`--run_flair false`, AEGIS and Mikado still receive the merged
Minimap2/StringTie long-read GTFs; the FLAIR input is just an empty sentinel.

## EDTA Repeat Annotation

Process: `EDTA`  
Module: `modules/EDTA.nf`  
Wrapper: `scripts/edta.sh`

```bash
bash scripts/edta.sh \
  -g <genome_fasta> \
  -n <cpus>
```

The wrapper runs:

```bash
export PYTHONNOUSERSITE=1

EDTA.pl \
  --genome <genome_fasta> \
  --species others \
  --step all \
  --sensitive 1 \
  --anno 1 \
  --overwrite 1 \
  --threads <cpus> \
  --force 1
```

Important options:

| Option | Purpose | Why it matters |
|---|---|---|
| `--species others` | Generic species mode. | Keeps the wrapper broadly usable across eukaryotes. |
| `--step all` | Run the complete EDTA workflow. | Produces library, annotation and masking outputs. |
| `--sensitive 1` | Sensitive repeat discovery. | Improves TE detection in complex genomes. |
| `--anno 1` | Annotate discovered repeats. | Produces the TE annotation GFF3. |
| `--overwrite 1` | Replace previous EDTA outputs in the task directory. | Keeps resumed task outputs deterministic. |
| `--threads` | Parallel execution. | Uses the task CPU allocation. |
| `--force 1` | Force execution despite existing intermediate state. | Avoids stale partial outputs. |

TITAN requires exactly one `*TElib.fa`, one `*TEanno.gff3` and one
`*MAKER.masked`, copying them to `edta.TElib.fa`, `edta.TEanno.gff3` and
`edta.MAKER.masked`.

## Protein FASTA Normalization

Process: `normalize_protein_fastas`  
Module: `modules/normalize_protein_fastas.nf`

```bash
python3 scripts/clean_protein_fasta_for_BRAKER3.py \
  <protein_fasta> \
  protein_<safe_organism>.fa \
  <safe_organism>
```

This internal cleaner prepares protein FASTA headers/sequences for BRAKER3.
The normalized FASTAs are collected and passed as `--prot_seq`.

## BRAKER3

Processes:

* `braker3_prediction`
* `braker3_prediction_with_long_reads`

Modules:

* `modules/braker3_prediction.nf`
* `modules/braker3_prediction_with_long_reads.nf`

Both call `scripts/run_braker3_prediction.sh`. The short-read branch collects
STAR BAMs only; the long-read branch also collects Minimap2 BAMs.

```bash
/BRAKER-3.0.8/scripts/braker.pl \
  --genome=<genome> \
  --bam=<comma_separated_short_and_optional_long_bams> \
  --prot_seq=<comma_separated_cleaned_protein_fastas> \
  --threads=<cpus> \
  --workingdir=<task_workdir> \
  --softmasking \
  --gff3 \
  --PROTHINT_PATH=/ProtHint-2.6.0/bin/ \
  --GENEMARK_PATH=/GeneMark-ETP \
  --AUGUSTUS_CONFIG_PATH=<task_workdir>/augustus_config \
  --AUGUSTUS_BIN_PATH=/Augustus/bin \
  --AUGUSTUS_SCRIPTS_PATH=/Augustus/scripts \
  --TSEBRA_PATH=/TSEBRA/bin
```

Important options:

| Option | Purpose | Why it matters |
|---|---|---|
| `--genome` | Target genome FASTA. | Prediction coordinate reference. |
| `--bam` | RNA-seq evidence BAMs. | Supplies expression and splice hints. |
| `--prot_seq` | Protein homology evidence. | Enables protein-supported prediction. |
| `--softmasking` | Treat lowercase sequence as masked. | Uses repeat masking appropriately. |
| `--gff3` | GFF3 output mode. | Feeds AEGIS and Mikado directly. |
| `--PROTHINT_PATH`, `--GENEMARK_PATH`, `--AUGUSTUS_*`, `--TSEBRA_PATH` | Explicit tool paths. | Avoids relying on ambiguous container PATH state. |

TITAN copies `Augustus/augustus.hints.gff3`, `GeneMark-ETP/genemark.gtf`,
`GeneMark-ETP/genemark_supported.gtf` and `braker.gff3`.

## Helixer

Process: `helixer_prediction`  
Module: `modules/helixer_prediction.nf`

Runs when `--run_helixer true`; otherwise TITAN emits an empty `helixer.gff3`.

```bash
export OMP_NUM_THREADS=<cpus>
export TF_NUM_INTRAOP_THREADS=<cpus>
export TF_NUM_INTEROP_THREADS=1

Helixer.py \
  --fasta-path <edta_masked_genome> \
  --lineage <params.helixer_model> \
  --gff-output-path helixer.gff3 \
  --temporary-dir helixer_tmp \
  --downloaded-model-path <params.helixer_model_dir>
```

Important options:

| Option | Purpose | Why it matters |
|---|---|---|
| `--lineage` | Selects the model lineage. | Must match the staged Helixer model. |
| `--downloaded-model-path` | Offline model location. | Keeps production runs reproducible. |
| TensorFlow thread env vars | Bound CPU threading. | Prevents uncontrolled CPU oversubscription. |

## AEGIS Evidence Merge

Process: `aegis_merge`  
Module: `modules/aegis_merge.nf`  
Wrapper: `scripts/run_aegis_merge.sh`

AEGIS is the central TITAN integration step. The wrapper first validates
required and optional evidence, then records included inputs in
`aegis_inputs.tsv`.

Required evidence in all modes:

* Liftoff GFF3
* BRAKER3 AUGUSTUS GFF3
* BRAKER3 GeneMark GTF
* EGAPx GFF3
* STAR/StringTie stranded default GTF
* STAR/StringTie stranded alternative GTF
* STAR/PsiCLASS stranded GTF

Required only when long reads are detected:

* Minimap2/StringTie long-read default GTF
* Minimap2/StringTie long-read alternative GTF

Optional evidence:

* STAR/PsiCLASS unstranded GTF
* STAR/StringTie unstranded default and alternative GTFs
* FLAIR isoforms GTF
* Helixer GFF3

Merge:

```bash
/opt/conda/envs/bio_env/bin/python -m aegis merge \
  -d aegis_merge \
  -o final_annotation \
  <included_evidence_files...>
```

Rename:

```bash
/opt/conda/envs/bio_env/bin/python -m aegis rename \
  -a final_annotation \
  -d aegis_rename \
  --prefix <params.aegis_gene_id_prefix> \
  --gene-id-correspondences \
  aegis_merge/final_annotation.gff3
```

Tidy:

```bash
/opt/conda/envs/bio_env/bin/python -m aegis tidy \
  -a final_annotation_renamed \
  -d aegis_tidy \
  --standard-features \
  aegis_rename/final_annotation_renamed.gff3
```

Extract proteins:

```bash
/opt/conda/envs/bio_env/bin/python -m aegis extract \
  -f protein \
  -m all \
  -d aegis_proteins_all \
  final_annotation.gff3 \
  <edta_masked_genome>

/opt/conda/envs/bio_env/bin/python -m aegis extract \
  -f protein \
  -m main \
  -d aegis_proteins_main \
  final_annotation.gff3 \
  <edta_masked_genome>
```

Important options:

| Option | Purpose | Why it matters |
|---|---|---|
| `merge -d` | AEGIS merge work/output directory. | Keeps merge products isolated. |
| `merge -o` | Annotation prefix. | Standardizes downstream file names. |
| `rename --prefix` | Final gene ID prefix. | Produces project-specific stable identifiers. |
| `rename --gene-id-correspondences` | ID mapping table. | Preserves traceability after renaming. |
| `tidy --standard-features` | Normalize feature structure. | Improves final GFF3 consistency. |
| `extract -m all` / `-m main` | Protein extraction mode. | Produces all and main protein sets for functional annotation and QC. |

## Diamond2GO

Process: `diamond2go`  
Module: `modules/diamond2go.nf`

Runs on both AEGIS protein FASTAs: all proteins and main proteins. TITAN wraps
the `diamond` binary so Diamond2GO uses the Nextflow CPU allocation.

```bash
perl /Diamond2GO/Diamond2go.pl \
  -d /Diamond2GO/resources/nr_clean_d2go.dmnd \
  -q query.fasta \
  -t protein
```

Important options:

| Option | Purpose | Why it matters |
|---|---|---|
| `-d` | Diamond2GO DIAMOND database. | Functional annotation search target. |
| `-q` | Query protein FASTA. | Runs separately for all and main proteins. |
| `-t protein` | Protein query type. | Matches AEGIS extracted proteins. |

## eggNOG-mapper

Process: `eggnog_mapper`  
Module: `modules/eggnog_mapper.nf`

Runs when `--run_eggnog_mapper true`; otherwise TITAN emits empty annotation
files.

```bash
emapper.py \
  -i <proteins.fa> \
  --itype proteins \
  -m diamond \
  --cpu <cpus> \
  --data_dir <params.eggnog_data_dir> \
  --output <final_annotation_proteins_all_or_main> \
  --output_dir . \
  --excel \
  --report_orthologs \
  --sensmode <params.eggnog_mapper_sensmode> \
  --override \
  --tax_scope <params.eggnog_mapper_tax_scope>
```

`--tax_scope` is added only when `--eggnog_mapper_tax_scope` is set.

Important options:

| Option | Purpose | Why it matters |
|---|---|---|
| `--itype proteins` | Protein input mode. | Matches final AEGIS protein FASTAs. |
| `-m diamond` | DIAMOND search backend. | Fast protein homology search. |
| `--data_dir` | Offline eggNOG data. | Required for reproducible production runs. |
| `--excel` | XLSX output. | User-facing spreadsheet output. |
| `--report_orthologs` | Ortholog reports. | Adds orthology detail beyond annotations. |
| `--sensmode` | Search sensitivity. | Tunes speed/sensitivity tradeoff. |
| `--override` | Replace existing output names. | Keeps resumed task directories deterministic. |

## InterProScan

Process: `interproscan`  
Module: `modules/interproscan.nf`

Runs when `--run_interproscan true`; otherwise TITAN emits empty TSV/GFF3/JSON
files.

```bash
/opt/interproscan/interproscan.sh \
  -i <proteins.fa> \
  -b <final_annotation_proteins_all_or_main> \
  -T interproscan_tmp \
  -cpu <cpus> \
  -dp \
  -f TSV,GFF3,JSON \
  -goterms \
  -pathways
```

Important options:

| Option | Purpose | Why it matters |
|---|---|---|
| `-b` | Output basename. | Separates all and main protein results. |
| `-T` | Temporary directory. | Keeps InterProScan scratch inside the task. |
| `-dp` | Disable precalculated match lookup. | Supports offline production runs. |
| `-f TSV,GFF3,JSON` | Output formats. | Provides tabular, feature and structured results. |
| `-goterms`, `-pathways` | Add GO and pathway annotations. | Enriches functional annotation outputs. |

## Mikado And TransDecoder

Process: `mikado_prepare`  
Module: `modules/mikado.nf`

Runs when `--run_mikado true`; otherwise TITAN emits empty/skipped sentinel
outputs.

```bash
awk -F'\t' 'NF == 9 || /^#/' <liftoff_gff3> > liftoff_sanitized.gff3

python3 scripts/make_mikado_list.py \
  --source liftoff_sanitized.gff3:liftoff:False:10:True \
  --source <egapx_gff3>:egapx:False:9:True \
  --source <braker_augustus_gff3>:braker_augustus:False:8:False \
  --source <braker_genemark_gtf>:braker_genemark:False:7:False \
  --source <star_stringtie_default_stranded>:star_stringtie_default_stranded:True:5:False \
  --source <star_stringtie_alt_stranded>:star_stringtie_alt_stranded:True:4:False \
  --source <star_psiclass_stranded>:star_psiclass_stranded:True:5:False \
  --source <star_psiclass_unstranded>:star_psiclass_unstranded:False:3:False \
  --source <star_stringtie_default_unstranded>:star_stringtie_default_unstranded:False:3:False \
  --source <star_stringtie_alt_unstranded>:star_stringtie_alt_unstranded:False:2:False \
  --source <hisat2_stringtie_default_stranded>:hisat2_stringtie_default_stranded:True:5:False \
  --source <hisat2_stringtie_alt_stranded>:hisat2_stringtie_alt_stranded:True:4:False \
  --source <hisat2_stringtie_default_unstranded>:hisat2_stringtie_default_unstranded:False:3:False \
  --source <hisat2_stringtie_alt_unstranded>:hisat2_stringtie_alt_unstranded:False:2:False \
  --source <long_reads_default>:long_reads_default:False:6:False \
  --source <long_reads_alt>:long_reads_alt:False:5:False \
  --source <flair_isoforms_gtf>:flair_isoforms:False:6:False \
  --source <helixer_gff3>:helixer:False:4:False \
  -o transcript_inputs.tsv

mikado configure \
  --list transcript_inputs.tsv \
  --reference <genome> \
  --mode <params.mikado_mode> \
  --scoring <params.mikado_scoring> \
  mikado_configuration.yaml

mikado prepare --json-conf mikado_configuration.yaml
```

`--source` entries use `path:label:stranded:score:is_reference`.

Process: `transdecoder_longorfs`  
Module: `modules/transdecoder.nf`

Runs when both `--run_mikado true` and `--run_transdecoder true`.

```bash
export PATH=/usr/local/opt/transdecoder/util:/usr/local/opt/transdecoder:$PATH
cp <mikado_prepared.fasta> mikado_prepared.fasta
TransDecoder.LongOrfs -t mikado_prepared.fasta
```

Process: `transdecoder_predict`

```bash
export PATH=/usr/local/opt/transdecoder/util:/usr/local/opt/transdecoder:$PATH
cp <mikado_prepared.fasta> mikado_prepared.fasta
cp -r <longorfs_dir> mikado_prepared.fasta.transdecoder_dir
TransDecoder.Predict \
  -t mikado_prepared.fasta \
  --single_best_only
```

Process: `mikado_serialise`

```bash
mikado serialise \
  --json-conf <mikado_configuration.yaml> \
  --orfs <transdecoder.bed> \
  --procs <cpus>
```

Process: `mikado_pick`

```bash
mikado pick \
  --json-conf <mikado_configuration.yaml> \
  --subloci-out mikado.subloci.gff3 \
  --loci-out mikado.loci.gff3 \
  --procs <cpus>
```

Process: `final_annotation_sources_qc`

```bash
python3 scripts/compare_final_annotations.py \
  --aegis-gff3 <aegis_final_annotation.gff3> \
  --mikado-gff3 <final_mikado_annotation.gff3> \
  --json-report final_annotation_sources.json \
  --multiqc-tsv final_annotation_sources_mqc.tsv
```

## lncRNA Candidate Annotation

Process: `lncrna_candidate_annotation`  
Module: `modules/lncrna_candidate_annotation.nf`

Runs when `--run_lncrna true`; otherwise TITAN emits empty candidate and CPAT
outputs.

Initial candidate build:

```bash
python3 scripts/build_lncrna_candidates.py \
  --genome <genome> \
  --final-annotation <final_annotation.gff3> \
  --trna-gff3 <trna.gff3> \
  --rfam-gff3 <rfam_ncrna.gff3> \
  --min-length <params.lncrna_min_length> \
  --output-prefix lncrna_candidates \
  <star_stringtie_gtf> <hisat2_stringtie_gtf> <long_reads_gtf>
```

Optional CPAT plant scoring when model files exist and candidates are present:

```bash
cpat.py \
  -x <params.cpat_model_dir>/Plant_Hexamer.tsv \
  -d <params.cpat_model_dir>/Plant.logit.RData \
  -g lncrna_candidates.fasta \
  -o cpat_plant.output
```

Final candidate rebuild with CPAT filtering:

```bash
python3 scripts/build_lncrna_candidates.py \
  --genome <genome> \
  --final-annotation <final_annotation.gff3> \
  --trna-gff3 <trna.gff3> \
  --rfam-gff3 <rfam_ncrna.gff3> \
  --min-length <params.lncrna_min_length> \
  --cpat-best-tsv cpat_plant.output.ORF_prob.best.tsv \
  --cpat-cutoff <params.cpat_plant_cutoff> \
  --output-prefix lncrna_candidates \
  <star_stringtie_gtf> <hisat2_stringtie_gtf> <long_reads_gtf>
```

Important options:

| Option | Purpose | Why it matters |
|---|---|---|
| `--final-annotation` | Final AEGIS coding annotation. | Excludes coding overlaps. |
| `--trna-gff3`, `--rfam-gff3` | ncRNA exclusion tracks. | Avoids labeling known ncRNAs as lncRNA candidates. |
| `--min-length` | Minimum candidate length. | Removes short fragments. |
| `--cpat-cutoff` | Coding probability threshold. | Filters likely coding transcripts. |

## SQANTI3 Long-Read QC

Processes:

* `sqanti3_qc`
* `sqanti3_qc_multiqc`

Module: `modules/sqanti3_qc.nf`

Runs when `--run_sqanti3 true` and long reads are present. TITAN runs it
separately for Minimap2/StringTie long-read transcripts and FLAIR isoforms.

```bash
gffread <final_annotation.gff3> \
  -T \
  -o reference.sqanti3.gtf

sqanti3_qc.py \
  <isoforms_gtf> \
  reference.sqanti3.gtf \
  <genome> \
  --dir sqanti3_out \
  --output <source_label> \
  --cpus <cpus>
```

Then TITAN summarizes structural categories from
`<source_label>.sqanti3_classification.txt`. The MultiQC process merges the two
summary TSV files into `sqanti3_long_read_mqc.tsv`.

Important options:

| Option | Purpose | Why it matters |
|---|---|---|
| `gffread -T` | Convert GFF3 reference to GTF. | SQANTI3 expects a reference GTF. |
| `--dir` | SQANTI3 output directory. | Keeps raw SQANTI3 products grouped. |
| `--output` | Source-specific basename. | Separates StringTie and FLAIR QC. |
| `--cpus` | Parallelism. | Matches task CPU allocation. |

## Final Annotation Validation

Process: `validate_final_annotation`  
Module: `modules/validate_final_annotation.nf`

```bash
python3 scripts/validate_final_annotation.py \
  --genome <masked_genome> \
  --annotation <final_annotation.gff3> \
  --proteins-all <final_annotation_proteins_all.fasta> \
  --proteins-main <final_annotation_proteins_main.fasta> \
  --json-report final_annotation_validation.json \
  --text-report final_annotation_validation.txt
```

This internal validator checks final GFF3/protein consistency before the run
is considered complete.

## Final Transcriptome And Expression Support

Process: `final_transcriptome_index`  
Module: `modules/final_expression_validation.nf`

Runs when `--run_expression_validation true`; otherwise TITAN emits skipped
sentinels.

```bash
python3 - <<'PY'
# Extract transcript sequences from <genome> and <final_annotation.gff3>
# into final_transcripts.fasta.
PY

salmon index \
  -t final_transcripts.fasta \
  -i final_salmon_index \
  -p <cpus>
```

Process: `final_expression_quant`

Paired-end:

```bash
salmon quant \
  -i <final_salmon_index> \
  -l A \
  -p <cpus> \
  -1 <read_1> \
  -2 <read_2> \
  -o <sample_ID>_quant \
  --validateMappings \
  2> <sample_ID>.expression_quant.log
```

Single-end:

```bash
salmon quant \
  -i <final_salmon_index> \
  -l A \
  -p <cpus> \
  -r <read_1> \
  -o <sample_ID>_quant \
  --validateMappings \
  2> <sample_ID>.expression_quant.log
```

Process: `expression_support_summary`

```bash
python3 scripts/summarize_expression_support.py \
  --gff <final_annotation.gff3> \
  --quant-dirs quant_dirs/* \
  --min-tpm <params.expression_support_min_tpm> \
  -o expression_support_summary.json \
  --multiqc-tsv expression_support_summary_mqc.tsv \
  --tpm-matrix gene_tpm_matrix.tsv
```

Important options:

| Option | Purpose | Why it matters |
|---|---|---|
| `-l A` | Autodetect library type. | Avoids forcing strandness during final support QC. |
| `--validateMappings` | Selective alignment validation. | Improves quantification reliability. |
| `--min-tpm` | Support threshold. | Defines supported versus unsupported genes. |

## BUSCO

Process: `busco`  
Module: `modules/busco.nf`

Runs when `--run_busco true`; otherwise TITAN emits skipped sentinel outputs.

```bash
busco \
  -i <final_annotation_proteins_main.fasta> \
  -m protein \
  -l <params.busco_lineage> \
  -o busco_out \
  -c <cpus> \
  --offline \
  --download_path <params.busco_data_dir> \
  -f
```

Important options:

| Option | Purpose | Why it matters |
|---|---|---|
| `-m protein` | Protein mode. | Evaluates the final AEGIS main protein set. |
| `-l` | BUSCO lineage. | Determines expected ortholog set. |
| `--offline`, `--download_path` | Use staged BUSCO data. | Required for reproducible production runs. |
| `-f` | Force output directory reuse. | Keeps task reruns deterministic. |

## OMArk

Process: `omark`  
Module: `modules/omark.nf`

Runs when `--run_omark true`; otherwise TITAN emits skipped sentinel outputs.

```bash
omamer search \
  --db <params.omark_data_dir>/omamer.h5 \
  --query <final_annotation_proteins_main.fasta> \
  --out proteins_main.omamer \
  --nthreads <cpus>

omark \
  -f proteins_main.omamer \
  -d <params.omark_data_dir>/omamer.h5 \
  -o omark_out
```

Important options:

| Option | Purpose | Why it matters |
|---|---|---|
| `omamer search --db` | OMAmer database. | Provides hierarchical orthology placement. |
| `--query` | Main protein FASTA. | QC is run on the representative final protein set. |
| `omark -f` | OMAmer search result. | Input for completeness/consistency checks. |
| `omark -d` | Same OMAmer database. | Keeps OMArk and OMAmer data synchronized. |

## AGAT Structural Statistics

Process: `agat_stats`  
Module: `modules/agat_stats.nf`

```bash
agat_sp_statistics.pl \
  --gff <final_annotation.gff3> \
  --output agat_stats.txt
```

This reports structural statistics for the final AEGIS GFF3.

## ncRNA QC

Process: `ncrna_annotation_qc`  
Module: `modules/ncrna_annotation_qc.nf`

```bash
agat_sp_statistics.pl \
  --gff <trna.gff3> \
  --output trna_agat_stats.txt

agat_sp_statistics.pl \
  --gff <rfam_ncrna.gff3> \
  --output rfam_agat_stats.txt

python3 - <<'PY'
# Count tRNA and Rfam feature records and write ncrna_mqc.tsv.
PY
```

This branch is QC/reporting only; tRNA and Rfam tracks are not merged into the
final AEGIS coding annotation.

## MultiQC

Process: `multiqc_report`  
Module: `modules/multiqc_report.nf`

```bash
mkdir -p mqc_input
cp <fastp_json_reports> mqc_input/
cp <busco_short_summary> mqc_input/
cp <omark_mqc_tsv> mqc_input/
cp <agat_stats_txt> mqc_input/
cp <ncrna_mqc_tsv> mqc_input/
cp <lncrna_mqc_tsv> mqc_input/
cp <sqanti3_mqc_tsv> mqc_input/
cp <expression_support_mqc_tsv> mqc_input/
cp <final_annotation_sources_mqc_tsv> mqc_input/
cp <final_annotation_validation_json> mqc_input/

multiqc \
  mqc_input \
  --filename titan_multiqc_report.html \
  --force
```

Important options:

| Option | Purpose | Why it matters |
|---|---|---|
| `mqc_input` | Staged report inputs. | Keeps MultiQC discovery deterministic. |
| `--filename` | Stable report name. | Predictable final HTML output. |
| `--force` | Overwrite existing report. | Supports resumed/re-run task directories. |

## Provenance Manifests

Processes:

* `additional_annotations_provenance`
* `titan_provenance`

Modules:

* `modules/additional_annotations_provenance.nf`
* `modules/titan_provenance.nf`

Both use inline Python to compute file records and SHA-256 checksums:

```bash
python3 - <<'PY'
# Build JSON manifest with paths, sizes, checksums, workflow metadata,
# parameters, output records and module version files.
PY
```

`additional_annotations_provenance` writes
`additional_annotations_manifest.json`. `titan_provenance` writes the global
`evidence_manifest.json`.

## Skipped And Sentinel Processes

Some processes deliberately emit empty placeholders instead of running a tool:

| Process | Trigger | Sentinel purpose |
|---|---|---|
| `empty_long_read_evidence` | No long-read samples. | Provides empty long-read GTFs for AEGIS/Mikado inputs. |
| `flair_empty_evidence` | No long reads or FLAIR disabled. | Provides empty FLAIR GTF/FASTA. |
| `trnascan_se` | `--run_trnascan false`. | Preserves ncRNA QC and lncRNA inputs. |
| `infernal_rfam_search` | `--run_rfam false`. | Preserves Rfam merge/QC inputs. |
| `helixer_prediction` | `--run_helixer false`. | Optional AEGIS/Mikado evidence is absent but typed. |
| `eggnog_mapper` | `--run_eggnog_mapper false`. | Functional output files remain predictable. |
| `interproscan` | `--run_interproscan false`. | Functional output files remain predictable. |
| `mikado_*`, `transdecoder_*` | Mikado/TransDecoder disabled. | Source-comparison branch remains stable. |
| `lncrna_candidate_annotation` | `--run_lncrna false`. | MultiQC and provenance still receive files. |
| `sqanti3_qc` | SQANTI3 disabled or no long reads. | MultiQC still receives zero-count summaries. |
| `final_expression_*` | `--run_expression_validation false`. | Expression support outputs remain predictable. |
| `busco`, `omark` | Their run flags are false. | QC report still records skipped status. |

## Process Coverage Inventory

This inventory lists every process declared in `modules/*.nf` and the section
that documents its command or sentinel behavior.

| Process | Module | Covered in section |
|---|---|---|
| `EDTA` | `modules/EDTA.nf` | EDTA Repeat Annotation |
| `Stringtie_merging_long_reads` | `modules/Stringtie_merging_long_reads.nf` | Merging Transcript Assemblies |
| `Stringtie_merging_short_reads_STAR` | `modules/Stringtie_merging_short_reads_STAR.nf` | Merging Transcript Assemblies |
| `Stringtie_merging_short_reads_hisat2` | `modules/Stringtie_merging_short_reads_hisat2.nf` | Merging Transcript Assemblies |
| `additional_annotations_provenance` | `modules/additional_annotations_provenance.nf` | Provenance Manifests |
| `aegis_merge` | `modules/aegis_merge.nf` | AEGIS Evidence Merge |
| `agat_convert_gff3_to_cds_fasta` | `modules/agat_convert_gff3_to_cds_fasta.nf` | AGAT Conversion For Liftoff Evidence |
| `agat_convert_gff3_to_gtf` | `modules/agat_convert_gff3_to_gtf.nf` | AGAT Conversion For Liftoff Evidence |
| `agat_stats` | `modules/agat_stats.nf` | AGAT Structural Statistics |
| `assembly_transcriptome_hisat2_stringtie` | `modules/assembly_transcriptome_hisat2_stringtie.nf` | StringTie Transcript Assembly |
| `assembly_transcriptome_minimap2_stringtie` | `modules/assembly_transcriptome_minimap2_stringtie.nf` | StringTie Transcript Assembly |
| `assembly_transcriptome_star_psiclass` | `modules/assembly_transcriptome_star_psiclass.nf` | PsiCLASS Transcript Assembly |
| `assembly_transcriptome_star_stringtie` | `modules/assembly_transcriptome_star_stringtie.nf` | StringTie Transcript Assembly |
| `braker3_prediction` | `modules/braker3_prediction.nf` | BRAKER3 |
| `braker3_prediction_with_long_reads` | `modules/braker3_prediction_with_long_reads.nf` | BRAKER3 |
| `busco` | `modules/busco.nf` | BUSCO |
| `clean_liftoff_gff3_for_agat` | `modules/clean_liftoff_gff3_for_agat.nf` | AGAT Conversion For Liftoff Evidence |
| `diamond2go` | `modules/diamond2go.nf` | Diamond2GO |
| `egapx` | `modules/egapx.nf` | EGAPx |
| `eggnog_mapper` | `modules/eggnog_mapper.nf` | eggNOG-mapper |
| `empty_long_read_evidence` | `modules/empty_long_read_evidence.nf` | Skipped And Sentinel Processes |
| `expression_support_summary` | `modules/final_expression_validation.nf` | Final Transcriptome And Expression Support |
| `final_annotation_sources_qc` | `modules/mikado.nf` | Mikado And TransDecoder |
| `final_expression_quant` | `modules/final_expression_validation.nf` | Final Transcriptome And Expression Support |
| `final_transcriptome_index` | `modules/final_expression_validation.nf` | Final Transcriptome And Expression Support |
| `flair_empty_evidence` | `modules/flair.nf` | Skipped And Sentinel Processes |
| `flair_isoforms` | `modules/flair.nf` | FLAIR Isoforms |
| `flair_merge_isoforms` | `modules/flair.nf` | FLAIR Isoforms |
| `gffcompare` | `modules/gffcompare.nf` | GffCompare For PsiCLASS |
| `helixer_prediction` | `modules/helixer_prediction.nf` | Helixer |
| `hisat2_alignment` | `modules/hisat2_alignment.nf` | HISAT2 Index And Alignment |
| `hisat2_genome_indices` | `modules/hisat2_genome_indices.nf` | HISAT2 Index And Alignment |
| `infernal_rfam_merge` | `modules/infernal_rfam.nf` | Rfam Annotation |
| `infernal_rfam_search` | `modules/infernal_rfam.nf` | Rfam Annotation |
| `interproscan` | `modules/interproscan.nf` | InterProScan |
| `liftoff_annotations` | `modules/liftoff_annotations.nf` | Liftoff Transfer |
| `lncrna_candidate_annotation` | `modules/lncrna_candidate_annotation.nf` | lncRNA Candidate Annotation |
| `mikado_pick` | `modules/mikado.nf` | Mikado And TransDecoder |
| `mikado_prepare` | `modules/mikado.nf` | Mikado And TransDecoder |
| `mikado_serialise` | `modules/mikado.nf` | Mikado And TransDecoder |
| `minimap2_alignment` | `modules/minimap2_alignment.nf` | Long-Read Staging And Minimap2 |
| `minimap2_genome_indices` | `modules/minimap2_genome_indices.nf` | Long-Read Staging And Minimap2 |
| `multiqc_report` | `modules/multiqc_report.nf` | MultiQC |
| `ncrna_annotation_qc` | `modules/ncrna_annotation_qc.nf` | ncRNA QC |
| `normalize_protein_fastas` | `modules/normalize_protein_fastas.nf` | Protein FASTA Normalization |
| `omark` | `modules/omark.nf` | OMArk |
| `prepare_RNAseq_fastq_files_long` | `modules/prepare_RNAseq_fastq_files_long.nf` | Long-Read Staging And Minimap2 |
| `prepare_RNAseq_fastq_files_short` | `modules/prepare_RNAseq_fastq_files_short.nf` | Short-Read Staging |
| `rfam_split_genome` | `modules/infernal_rfam.nf` | Rfam Annotation |
| `salmon_index` | `modules/salmon_index.nf` | Salmon Strand Inference |
| `salmon_strand_inference` | `modules/salmon_strand_inference.nf` | Salmon Strand Inference |
| `sqanti3_qc` | `modules/sqanti3_qc.nf` | SQANTI3 Long-Read QC |
| `sqanti3_qc_multiqc` | `modules/sqanti3_qc.nf` | SQANTI3 Long-Read QC |
| `star_alignment` | `modules/star_alignment.nf` | STAR Index And Alignment |
| `star_genome_indices` | `modules/star_genome_indices.nf` | STAR Index And Alignment |
| `titan_provenance` | `modules/titan_provenance.nf` | Provenance Manifests |
| `transdecoder_longorfs` | `modules/transdecoder.nf` | Mikado And TransDecoder |
| `transdecoder_predict` | `modules/transdecoder.nf` | Mikado And TransDecoder |
| `trimming_fastq` | `modules/trimming_fastq.nf` | Short-Read Staging |
| `trnascan_se` | `modules/trnascan_se.nf` | tRNA Annotation |
| `trnascan_to_gff3` | `modules/trnascan_se.nf` | tRNA Annotation |
| `validate_final_annotation` | `modules/validate_final_annotation.nf` | Final Annotation Validation |
| `validate_inputs` | `modules/validate_inputs.nf` | Input Validation |
