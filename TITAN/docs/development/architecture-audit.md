# TITAN architecture audit

Date: 2026-07-16
Scope: static audit of `main.nf`, `subworkflows/*.nf`, `modules/*.nf`, `conf/*.config`, helper scripts and current development docs.

This audit focuses on making TITAN a maintainable Nextflow workflow before adding more scientific surface area. It intentionally does not validate biological results.

## Executive summary

TITAN has a useful biological decomposition, but the workflow is not yet architecturally solid. The main problems are not isolated syntax issues; they are contract issues between workflow layers.

Highest-risk findings:

* `main.nf` still acts as orchestration, validation and Aegis evidence loader instead of being a thin entrypoint.
* `generate_evidence_data` emits a string-keyed mixed channel, while `aegis` expects a manually reconstructed list of files. This bypasses normal Nextflow dataflow guarantees.
* `aegis` mode reads files back from `params.output_dir` by filename. That means Aegis is coupled to publish side effects instead of process outputs.
* EDTA is now mandatory in the evidence workflow, and Aegis-only runs fail early if the hard-masked genome is missing. The remaining issue is that Aegis still discovers that file from `output_dir` instead of a typed channel or manifest.
* EGAPx is now mandatory in `generate_evidence_data` and its module input is staged as `path egapx_paramfile`. Its outputs are still captured broadly and are not yet named or consumed by Aegis.
* Many modules take `val()` inputs and then ignore them, reading from mounted `params.output_dir`, `projectDir/data/*` or helper scripts instead. This makes caching, resume and portability fragile.
* Long-read processing is now detected from `RNAseq_samplesheet` rows where `library_layout` is `long`; the old `use_long_reads` flag is no longer part of the runtime contract.
* Multiple modules copy files manually to `/outputdir` in addition to `publishDir`, creating hidden dependencies and duplicate output semantics.
* Several processes use `containerOptions` mounts for Docker-specific paths. This fights Nextflow staging and weakens Apptainer/HPC portability.

The recommended path is to stabilize contracts before changing scientific behavior:

1. Create a thin `main.nf` and a canonical `workflows/titan.nf`.
2. Introduce typed evidence outputs rather than string-key/file lists.
3. Make Aegis consume process outputs directly, not files discovered in `output_dir`.
4. Move the required EDTA hard-masked genome from `output_dir` discovery to a typed Aegis input or evidence manifest.
5. Replace broad EGAPx output capture with named EGAPx emits and add its GFF3 to the same evidence contract once output names are confirmed.
6. Migrate modules to staged `path` inputs and named emits module by module.

## Current topology

Current control flow:

```text
main.nf
  validates parameters
  builds RNA/protein channels for generate_evidence_data
  OR scans output_dir for Aegis evidence files

subworkflows/generate_evidence_data.nf
  prepares reads
  runs Liftoff, AGAT, Salmon, STAR, HISAT2, Minimap2, StringTie, PsiCLASS, GFFCompare, EDTA, BRAKER3
  emits one mixed channel of [string_key, file] pairs

subworkflows/aegis.nf
  receives a list of [string_key, file] pairs
  reconstructs arrays by key
  conditionally runs aegis_short_reads or aegis_long_reads
  conditionally runs Diamond2GO
```

Target control flow:

```text
main.nf
  include workflow TITAN
  no biological orchestration

workflows/titan.nf
  validate params
  build typed input channels
  run evidence generation if requested
  run Aegis integration if requested

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

### 1. Entrypoint and orchestration are mixed

Evidence:

* `main.nf` includes subworkflows, mutates defaults, validates parameters, builds CSV channels, scans `output_dir`, creates placeholder files and dispatches Aegis.
* `params.workflow == "all"` was accepted by validation but errored later before P1-001; it is now rejected early until a real `all` mode exists.

Impact:

* hard to test one layer at a time;
* hard to add EGAPx cleanly;
* hard to support a true `all` mode;
* user-facing workflow modes are implementation details rather than stable workflow contracts.

Recommendation:

* keep `main.nf` as a thin entrypoint;
* move orchestration into `workflows/titan.nf`;
* split Aegis-only and full-run paths at workflow level, not by scanning output directories.

### 2. Aegis is coupled to published files instead of process outputs

Evidence:

* `main.nf` in `--workflow aegis` checks hard-coded names such as `assembly_masked.EDTA.fasta`, `augustus.hints.gff3`, `genemark.gtf`, `liftoff_previous_annotations.gff3`, `merged_star_stringtie_stranded_default.gtf` inside `params.output_dir`.
* Missing required files are printed as `ERROR`, but execution continues until Aegis is skipped or receives incomplete arrays.

Impact:

* `-resume` and caching semantics are weakened;
* a file copied manually by a script becomes part of the workflow contract;
* a failed evidence process can be hidden if an old file exists in `output_dir`;
* Aegis cannot be reliably unit-tested as a subworkflow.

Recommendation:

* define an evidence bundle contract as named channels;
* make Aegis consume those channels directly for full runs;
* for Aegis-only reruns, introduce an explicit `--evidence_manifest` rather than ad hoc filename discovery.

### 3. EDTA is mandatory but still coupled to published filenames

Evidence:

* `generate_evidence_data` now always runs EDTA.
* `aegis` now always requires `masked_genome.masked_genome` and fails if it is missing.
* In Aegis-only mode, `main.nf` still discovers `assembly_masked.EDTA.fasta` by filename in `params.output_dir`.
* `EDTA.nf` emits `*MAKER.masked` but also manually copies it to `/outputdir/assembly_masked.EDTA.fasta`.

Impact:

* Aegis requirements are now explicit but still file-name based in Aegis-only mode;
* EDTA output naming depends on manual copy side effects.

Recommendation:

* Aegis requires exactly one hard-masked genome input;
* replace Aegis-only file discovery with `--evidence_manifest` or explicit `--masked_assembly`;
* emit `masked_genome` with the final public name directly from the EDTA process.

### 4. EGAPx is integrated in evidence generation but not typed

Evidence:

* `generate_evidence_data` now includes and calls `egapx(file(params.egapx_paramfile))`.
* `modules/egapx.nf` declares `path egapx_paramfile` and emits `egapx_results`.
* output is still `path("*")`, unnamed and too broad for stable Aegis integration.

Impact:

* EGAPx is mandatory in evidence generation;
* the output cannot yet be connected cleanly to Aegis;
* broad output capture can publish unrelated files and makes tests brittle.

Recommendation:

* emit named outputs such as `egapx_gff3`, `egapx_proteins`, `egapx_report` once actual EGAPx output names are confirmed;
* integrate `egapx_gff3` into the evidence bundle and Aegis command line when the output contract is stable.

### 5. Module inputs are not consistently staged

Evidence:

* read preparation and trimming modules mount `${projectDir}/data/RNAseq_data` and reference `/RNAseq_data/${sample_ID}*.fastq.gz`.
* BRAKER3 and Aegis mount `${projectDir}/data/protein_data`, even though the current samplesheet can point to arbitrary fixture paths.
* merge modules read from mounted output directories with helper scripts instead of consuming the GTF files passed in their input value.
* `braker3_prediction` takes `val(bam_short)`, but then retrieves BAMs from `/alignments` by scanning directories.

Impact:

* Nextflow cannot track the real file dependencies;
* modules are difficult to test with minimal fixtures;
* moving to Apptainer or a different working directory is risky;
* stale files in published output directories can affect new runs.

Recommendation:

* replace directory scans with staged `path` inputs;
* carry sample metadata as maps or tuples;
* make every process command use its declared inputs;
* use `publishDir` only for public outputs, not as an internal data bus.

### 6. Channel shapes and optional branches need formal contracts

Evidence:

* `generate_evidence_data` emits a mixed channel of key/file pairs.
* optional unstranded outputs are emitted as optional paths, then mixed unconditionally.
* `aegis.nf` converts a list of pairs into arrays and indexes `[0]`.
* Before P1-001, long-read branch checks were split between `main.nf`, subworkflow `if (params.use_long_reads)`, and process-level `when: params.use_long_reads`. The branch is now controlled by samplesheet detection: at least one `library_layout=long` row enables long-read modules.

Impact:

* missing evidence can become an index error or a skipped process later;
* optional outputs are hard to reason about;
* boolean behavior can differ by layer;
* true `all` mode will be fragile until contracts are explicit.

Recommendation:

* define typed emits per evidence family;
* normalize booleans once in params validation;
* pass normalized booleans or channels to subworkflows;
* create explicit placeholder channels only for genuinely optional evidence, not required evidence.

### 7. Outputs and provenance are not cleanly separated

Evidence:

* EDTA, StringTie merging and GFFCompare copy public files manually to `/outputdir`.
* `publishDir` also publishes process outputs.
* output file names are sometimes command-dependent and sometimes manually renamed.
* no `versions.yml`, input manifest or evidence manifest is emitted.

Impact:

* public outputs are hard to document and validate;
* reruns can see old outputs;
* provenance cannot be audited automatically.

Recommendation:

* publish from process outputs only;
* emit an `evidence_manifest.json` for Aegis;
* emit `versions.yml` per module or per workflow;
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
3. Create `workflows/titan.nf` and keep `main.nf` thin.
4. Replace the mixed evidence channel with named emits.
5. Rewrite Aegis integration to consume named evidence directly.
6. Make EDTA hard-masked genome an explicit required Aegis input.
7. Replace broad mandatory EGAPx results with named EGAPx evidence emits and connect them to Aegis.
8. Migrate modules away from mounted output directories and hidden scans.
9. Add `-stub-run` and nf-test coverage for critical subworkflows.
10. Add CI once the contracts are stable.

## Immediate non-code acceptance criteria

Before large refactors, the following should be true:

* roadmap names the contract problems explicitly;
* one architecture document defines target channel shapes;
* every future implementation ticket has a rollback strategy;
* P0 test profile remains green after each step.
