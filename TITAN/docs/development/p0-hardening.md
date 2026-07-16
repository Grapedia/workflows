# TITAN P0 hardening summary

Date: 2026-07-16
Branch: `codex/titan-hardening`

This document summarizes the completed P0 work. It is a short operational guide; the detailed inventory and command log remain in `docs/development/audit.md`. The broader architecture audit and refactor target are in `docs/development/architecture-audit.md`.

## Scope

Completed items:

* `TITAN-P0-001`: establish a functional baseline.
* `TITAN-P0-002`: inventory inputs, outputs, tools and containers.
* `TITAN-P0-003`: add a minimal local test profile.
* `TITAN-P0-004`: create a minimal multi-input synthetic test dataset.
* `TITAN-P0-005`: validate required parameters before channel construction.

## Current baseline

The pipeline now resolves with Nextflow and the local test profile can run the full evidence-generation graph in stub mode without Slurm, Docker or production data:

```bash
nextflow config -profile test
nextflow run main.nf -profile test -stub-run -ansi-log false
```

The `test` profile uses:

* `conf/test.config`;
* local executor;
* Docker disabled;
* `test-data/minimal/valid` inputs;
* `test-results/` and `test-work/` transient directories;
* EDTA and EGAPx are mandatory in TITAN.
* Long reads are detected from `library_layout=long` in the RNA-seq samplesheet.

This is a workflow/bootstrap validation only. It does not validate the biological behavior of EDTA, EGAPx, BRAKER3, STAR, HISAT2, Minimap2, PsiCLASS, Liftoff or Aegis.

## Minimal test data

The synthetic fixture set is under `test-data/minimal`.

It covers:

* reference and target genome FASTA files;
* previous GFF3 annotation for Liftoff;
* RNA-seq samplesheet with `single`, `paired` and `long` layouts;
* tiny gzip-compressed FASTQ files matching the samplesheet;
* protein samplesheet with two protein evidence sources;
* EGAPx YAML parameter shape;
* precomputed minimal evidence files for EDTA, Liftoff, BRAKER/AUGUSTUS, GeneMark, STAR/StringTie, STAR/PsiCLASS and Minimap2/StringTie;
* negative fixtures for invalid FASTA/GFF3 validation.

Validate the fixture set with:

```bash
python3 scripts/validate_minimal_test_data.py
sha256sum -c test-data/minimal/checksums.sha256
```

## Parameter validation

`workflows/titan.nf` validates these required parameters inside the `TITAN` workflow before input channels are built:

```text
output_dir
egapx_paramfile
RNAseq_samplesheet
protein_samplesheet
new_assembly
previous_assembly
previous_annotations
```

The validation catches:

* missing or empty required values;
* flags passed without values, which Nextflow can expose as `true`;
* missing input files;
* the deprecated `--workflow` parameter, because TITAN now always runs the full graph.
The workflow also uses `Channel.fromPath(..., checkIfExists: true)` for samplesheets.

Negative checks used during P0-005:

```bash
nextflow run main.nf -profile test --RNAseq_samplesheet '' -stub-run -ansi-log false
nextflow run main.nf -profile test --previous_annotations test-data/minimal/valid/missing.gff3 -stub-run -ansi-log false
nextflow run main.nf -profile test --workflow aegis -stub-run -ansi-log false
```

Expected messages include:

```text
Missing required parameter(s): --RNAseq_samplesheet
Required input file(s) not found
--workflow is no longer supported
```

## Known limits after P0

The P0 work intentionally avoids large architectural refactors. Remaining issues for P1 and later:

* `workflows/titan.nf` still carries orchestration logic that should move progressively toward typed evidence contracts.
* Several images still use `latest`.
* EGAPx is mandatory in TITAN, emits named outputs, and `egapx_gff3` is consumed by AEGIS as an additional merge evidence.
* The test profile validates channel/file contracts in stub mode; it does not run scientific containers or heavy tools.
* There is no CI yet.
* Protein-related module volume mounts still assume production data layout in places such as `data/protein_data`.

P1-001 initially normalized the historical booleans; the current contract removes those biological switches and public partial workflow modes. EDTA and EGAPx are mandatory in TITAN, Aegis always runs after evidence generation, and long reads are inferred from the samplesheet. The remaining architectural work is tracked in `roadmap.md`.

## Relevant files

* `roadmap.md`: P0 status and acceptance criteria.
* `docs/development/audit.md`: detailed inventory, command log and risks.
* `conf/test.config`: local test profile.
* `test-data/minimal/README.md`: fixture-level documentation.
* `scripts/validate_minimal_test_data.py`: static fixture validator.
