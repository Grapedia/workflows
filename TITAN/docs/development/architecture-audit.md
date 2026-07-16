# TITAN architecture audit

Date: 2026-07-16
Scope: static audit of `main.nf`, `subworkflows/*.nf`, `modules/*.nf`, `conf/*.config`, helper scripts and current development docs.

This audit focuses on making TITAN a maintainable Nextflow workflow before adding more scientific surface area. It intentionally does not validate biological results.

## Executive summary

TITAN has a useful biological decomposition, but the workflow is not yet architecturally solid. The main problems are not isolated syntax issues; they are contract issues between workflow layers.

Highest-risk findings:

* P1-002 made `main.nf` a thin entrypoint; orchestration now lives in `workflows/titan.nf`.
* P1-003 replaced the string-keyed mixed evidence channel with named emits and explicit Aegis inputs.
* P1-004 removed public partial workflow modes: TITAN now always runs evidence generation and Aegis in one graph.
* Aegis consumes generated evidence directly through named channels, including the EDTA `masked_genome` output.
* EDTA is now mandatory in TITAN. The hard-masked genome is passed to Aegis as a typed channel.
* EGAPx is now mandatory, uses the official digest-pinned `ncbi/egapx` image through the official `v0.5.2` runner, emits named outputs, and its GFF3 is consumed by AEGIS as an additional merge evidence.
* P1-008 centralizes runtime containers in `nextflow.config` and locks them with `image@sha256` digests.
* P1-006 removed the highest-risk internal directory scans from evidence integration modules. StringTie merges, GFFCompare, BRAKER3, AGAT, Salmon and STAR/HISAT2/Minimap2 alignment stages now consume declared upstream files instead of rediscovering them from published directories.
* Long-read processing is now detected from `RNAseq_samplesheet` rows where `library_layout` is `long`; the old `use_long_reads` flag is no longer part of the runtime contract.
* P1-007 removed the remaining manual `/outputdir` publication from the core evidence path. Backward-compatible public names are now handled through declared outputs and `publishDir saveAs`.
* Several processes use `containerOptions` mounts for Docker-specific paths. This fights Nextflow staging and weakens Apptainer/HPC portability.

The recommended path is to stabilize contracts before changing scientific behavior:

1. Done in P1-002: keep `main.nf` thin and continue building the canonical orchestration in `workflows/titan.nf`.
2. Done in P1-003: introduce named evidence outputs and explicit Aegis inputs.
3. Done in P1-004: make Aegis consume process outputs directly, not files discovered in `output_dir`.
4. Done in P1-004: move the required EDTA hard-masked genome from `output_dir` discovery to a typed Aegis input.
5. Done in P1-005 and the follow-up AEGIS integration: replace broad EGAPx output capture with named EGAPx emits and pass `egapx_gff3` to `aegis merge`.
6. Done in P1-007 for provenance: publish `provenance/evidence_manifest.json` and `provenance/versions.yml`.

## Current topology

Current control flow:

```text
main.nf
  sets public fallback params
  includes and calls TITAN

workflows/titan.nf
  validates parameters
  builds RNA/protein channels
  runs evidence generation
  feeds named evidence channels into Aegis

subworkflows/generate_evidence_data.nf
  prepares reads
  runs Liftoff, AGAT, Salmon, STAR, HISAT2, Minimap2, StringTie, PsiCLASS, GFFCompare, EDTA, BRAKER3
  emits named evidence channels

subworkflows/aegis.nf
  receives explicit evidence inputs
  conditionally runs aegis_short_reads or aegis_long_reads
  conditionally runs Diamond2GO
```

Target control flow after remaining P1 work:

```text
main.nf
  include workflow TITAN
  no biological orchestration

workflows/titan.nf
  validate params
  build typed input channels
  run evidence generation
  run Aegis integration

subworkflows/evidence_generation
  emit named channels:
    liftoff_gff3
    egapx_gff3
    edta_masked_genome
    braker_augustus_gff3
    braker_genemark_gtf
    star_stringtie_stranded_default_gtf
    star_stringtie_stranded_alt_gtf
    star_psiclass_stranded_gtf
    optional unstranded evidence
    optional long-read evidence

subworkflows/aegis_integration
  take named evidence channels
  fail early if required evidence is missing
  run Aegis and Diamond2GO
```

## Findings

### 1. Entrypoint is thin; orchestration still needs contract cleanup

Evidence:

* `main.nf` now only sets public fallback params, includes `TITAN` from `workflows/titan.nf` and calls it.
* `workflows/titan.nf` includes subworkflows, validates parameters, builds CSV channels, dispatches `generate_evidence_data`, and passes its named outputs to Aegis.
* `params.workflow` is no longer part of the public API.

Impact:

* hard to test one layer at a time;
* hard to add EGAPx cleanly;
* partial reruns need an explicit future contract instead of hidden workflow modes.

Recommendation:

* keep `main.nf` thin;
* continue reducing orchestration complexity in `workflows/titan.nf`;
* keep the public workflow contract single-purpose; reintroduce partial reruns only with an explicit evidence manifest or resume strategy.

### 2. Aegis now consumes generated evidence directly

Evidence:

* `workflows/titan.nf` passes named emits from `generate_evidence_data` directly into `aegis`.
* The old Aegis-only path that checked hard-coded names in `params.output_dir` has been removed from the public workflow.

Impact:

* `-resume` and caching semantics are clearer for full TITAN runs;
* stale files in `output_dir` can no longer satisfy Aegis inputs accidentally;
* Aegis can be tested through its explicit input contract.

Recommendation:

* continue defining the evidence bundle as named channels;
* for future partial reruns, introduce an explicit `--evidence_manifest` rather than ad hoc filename discovery.

### 3. EDTA is mandatory and exposes its public masked-genome name declaratively

Evidence:

* TITAN now always runs EDTA.
* `aegis` now always requires `masked_genome.masked_genome` and fails if it is missing.
* `workflows/titan.nf` passes `evidence_data.masked_genome` directly to `aegis`.
* `EDTA.nf` emits `*MAKER.masked` and publishes the backward-compatible `assembly_masked.EDTA.fasta` name through `publishDir saveAs`.

Impact:

* Aegis requirements are explicit and channel-based;
* EDTA output naming is now declared in Nextflow publication rules.

Recommendation:

* Aegis requires exactly one hard-masked genome input;
* emit `masked_genome` with the final public name directly from the EDTA process.

### 4. EGAPx is integrated, typed and consumed by AEGIS

Evidence:

* `generate_evidence_data` now includes and calls `egapx(file(params.egapx_paramfile))` as an internal TITAN stage.
* `modules/egapx.nf` declares `path egapx_paramfile`.
* The module uses official EGAPx runner `v0.5.2` and the digest-pinned official `ncbi/egapx` Docker image.
* The module emits named EGAPx outputs: GFF3, GTF, proteins, CDS, transcripts, ASN, full output directory and versions.
* `workflows/titan.nf` passes `evidence_data.egapx_gff3` into `subworkflows/aegis.nf`.
* `aegis_short_reads` and `aegis_long_reads` include EGAPx GFF3 in the `aegis merge` input list.

Impact:

* EGAPx is mandatory in evidence generation;
* EGAPx output can now be referenced by stable channel names;
* AEGIS now treats EGAPx as an additional annotation evidence.

Recommendation:

* keep EGAPx as a named evidence input and document any future changes to merge priority explicitly.

### 5. Module inputs are increasingly staged, but cleanup is not complete

Evidence:

* `prepare_RNAseq_fastq_files_*` now materialize local/SRA reads in the workdir and emit staged read files.
* `trimming_fastq`, `minimap2_alignment`, `salmon_index`, `salmon_strand_inference`, `agat_convert_gff3_to_cds_fasta`, `star_alignment` and `hisat2_alignment` consume declared `path` inputs.
* StringTie merge modules and `gffcompare` receive explicit lists of staged GTF files from channel transformations.
* BRAKER3 receives staged genome, protein FASTA files and BAM lists instead of scanning `/alignments` or `${projectDir}/data/protein_data`.
* Local FASTQ/FASTA discovery is still convention-based through the explicit `RNAseq_data_dir` parameter and `sample_ID` naming.

Impact:

* Nextflow can now track the main evidence dependencies through staged inputs;
* stub tests exercise the same channel contracts used by AEGIS;
* stale files in published output directories no longer satisfy StringTie/GFFCompare/BRAKER3 inputs;
* RNA-seq local file naming still needs a future samplesheet schema upgrade if arbitrary file paths are required.

Recommendation:

* continue replacing remaining directory conventions with staged `path` inputs;
* carry sample metadata as maps or tuples;
* make every process command use its declared inputs;
* use `publishDir` only for public outputs, not as an internal data bus.

### 6. Channel shapes and optional branches need formal contracts

Evidence:

* P1-003 replaced the mixed key/file channel with named emits such as `masked_genome`, `liftoff_gff3`, `braker_augustus_gff3`, `star_stringtie_stranded_default_gtf` and long-read emits.
* optional unstranded outputs remain optional paths and should eventually get clearer placeholder/absence semantics.
* `aegis.nf` now takes explicit evidence inputs instead of converting a list of pairs into arrays.
* Before P1-001, long-read branch checks were split between `main.nf`, subworkflow `if (params.use_long_reads)`, and process-level `when: params.use_long_reads`. The branch is now controlled by samplesheet detection: at least one `library_layout=long` row enables long-read modules.

Impact:

* missing required evidence is represented by missing upstream process outputs;
* optional outputs are hard to reason about;
* boolean behavior can differ by layer;
* optional evidence still needs a cleaner absence model.

Recommendation:

* continue refining typed emits per evidence family;
* normalize booleans once in params validation;
* pass normalized booleans or channels to subworkflows;
* create explicit placeholder channels only for genuinely optional evidence, not required evidence.

### 7. Outputs and provenance are now explicitly published

Evidence:

* EDTA, StringTie merging and GFFCompare preserve historical public names with `publishDir saveAs`.
* `modules/titan_provenance.nf` emits `evidence_manifest.json` with input, evidence and final-output checksums.
* `modules/titan_provenance.nf` emits a workflow-level `versions.yml` for TITAN provenance, EGAPx and AEGIS versions.

Impact:

* public outputs are easier to document and validate;
* the final AEGIS evidence bundle can be audited after a run;
* module-level versions still need broader coverage in later work.

Recommendation:

* keep publishing from declared process outputs;
* extend `versions.yml` collection to remaining tools as container versions are locked;
* keep backward-compatible filenames through `publishDir saveAs` or stable process output names.

### 8. Profiles and resources are only partially normalized

Evidence:

* `conf/base.config` has generic labels, but modules hardcode many `cpus`.
* `nextflow.config` still has a global default of 100GB/20 CPU for `default`.
* profiles resolve, but production profiles still need validation with real Docker/Apptainer/Slurm behavior.

Impact:

* resource tuning requires code edits;
* local tests and HPC runs can diverge;
* some Docker `containerOptions` are not portable to Apptainer.

Recommendation:

* add labels per process class: indexing, alignment, transcript assembly, merge, prediction, annotation integration, lightweight validation;
* move resources into config;
* avoid Docker-specific mounts where Nextflow staging can handle inputs.

## Suggested migration order

The safest order is contract-first:

1. Normalize booleans and workflow mode semantics.
2. Add a samplesheet/schema validator for RNA-seq and proteins.
3. Done in P1-002: create `workflows/titan.nf` and keep `main.nf` thin.
4. Done in P1-003: replace the mixed evidence channel with named emits.
5. Done in P1-004: rewrite Aegis integration to consume named evidence directly.
6. Done in P1-004: make EDTA hard-masked genome an explicit required Aegis input.
7. Done in P1-005 and follow-up AEGIS integration: replace broad mandatory EGAPx results with named EGAPx evidence emits and connect `egapx_gff3` to AEGIS.
8. Done in P1-006 for the critical evidence path: migrate StringTie merges, GFFCompare, BRAKER3, AGAT, Salmon and STAR/HISAT2/Minimap2 alignment away from mounted output directories and hidden scans.
9. Add `-stub-run` and nf-test coverage for critical subworkflows.
10. Add CI once the contracts are stable.

## Immediate non-code acceptance criteria

Before large refactors, the following should be true:

* roadmap names the contract problems explicitly;
* one architecture document defines target channel shapes;
* every future implementation ticket has a rollback strategy;
* P0 test profile remains green after each step.
