# TITAN Architecture

## Table Of Contents

- [Workflow Shape](#workflow-shape)
- [Evidence Contracts](#evidence-contracts)
- [Mandatory Core Evidence](#mandatory-core-evidence)
- [RNA-Seq And Long Reads](#rna-seq-and-long-reads)
- [Additional Annotation Branches](#additional-annotation-branches)
- [Staging And Publication](#staging-and-publication)
- [Containers And Profiles](#containers-and-profiles)
- [Validation And Test Data](#validation-and-test-data)
- [Provenance](#provenance)


This document records durable architecture decisions for TITAN. It is written
for contributors who need to understand the workflow structure and the
contracts between layers.

## Workflow Shape

`main.nf` is intentionally a thin entrypoint. It defines public fallback
parameters, includes the `TITAN` workflow from `workflows/titan.nf`, and calls
it. Biological orchestration belongs in `workflows/titan.nf` and lower-level
subworkflows/modules.

The public runtime contract is a single full TITAN graph. TITAN no longer
exposes partial public workflow modes. Evidence generation, AEGIS integration,
additional annotation branches and reporting are connected in one graph so
Nextflow can track dependencies and cache behavior consistently.

Partial reruns should be introduced only through an explicit future contract,
such as an evidence manifest, not by rediscovering files from published output
directories.

## Evidence Contracts

Evidence is passed between workflow layers as named channels and typed process
outputs. Published output directories are for users, not an internal data bus.

Important evidence emits include:

* `masked_genome`
* `liftoff_gff3`
* `egapx_gff3`
* `braker_augustus_gff3`
* `braker_genemark_gtf`
* STAR/StringTie GTF evidence
* optional HISAT2/StringTie GTF evidence
* STAR/PsiCLASS GTF evidence
* optional long-read StringTie evidence
* optional FLAIR isoform evidence
* optional Helixer evidence

Mikado, not AEGIS, consumes these explicit evidence inputs from the workflow
graph and reconciles them into one gene set (see
[Additional Annotation Branches](#additional-annotation-branches)). AEGIS
only takes Mikado's consolidated output, `masked_genome` (for rename/tidy and
protein extraction) and `liftoff_gff3` + `previous_annotations` (for its
liftoff gene-ID carryover step). None of these layers use stale files
discovered under `params.output_dir`.

Optional branches should expose clear channel contracts. Required evidence
should fail early when absent. Optional evidence should use explicit absence or
sentinel semantics where downstream tools need stable inputs.

## Mandatory Core Evidence

EDTA is mandatory in TITAN. The EDTA hard-masked genome is a required channel
input to both Mikado (as the reference genome for `mikado configure`) and
AEGIS (for `aegis rename`/`tidy`/`extract`). The process publishes
backward-compatible filenames through `publishDir saveAs`, while downstream
workflow logic consumes declared outputs.

EGAPx is mandatory in evidence generation. TITAN runs the official EGAPx
runner with a digest-pinned EGAPx container, emits named EGAPx outputs, and
passes the EGAPx GFF3 to Mikado as annotation evidence (`source_score=9`,
`is_reference=True` in `mikado_prepare`).

Liftoff, transcriptome assemblies and BRAKER3 predictions form the rest of
the core evidence set consumed by Mikado. `run_mikado` and `run_transdecoder`
both default to `true` and are effectively required: `aegis_merge` fails fast
if `mikado_pick`'s output has no feature records. Liftoff's GFF3 is also
consumed a second time, directly by `aegis_liftoff_gene_ids`, to carry old
gene IDs over onto matching genes in AEGIS's output (see
[Additional Annotation Branches](#additional-annotation-branches)).

## RNA-Seq And Long Reads

RNA-seq samples are described by `RNAseq_samplesheet`. Short-read rows are
prepared, trimmed, aligned with STAR, assembled with StringTie and PsiCLASS,
and used for strand inference and final expression support. HISAT2/StringTie is
an opt-in evidence branch (`--run_hisat2 true`) that is disabled by default and
feeds Mikado, lncRNA candidates and provenance only; it does not feed AEGIS.

Long-read behavior is inferred from rows where `library_layout` is `long`.
There is no separate public `use_long_reads` switch. When long reads are
present, TITAN stages them, aligns them with Minimap2, assembles long-read
transcripts with StringTie, and can pass long-read evidence to Mikado, FLAIR
and SQANTI3.

## Additional Annotation Branches

tRNAscan-SE and Infernal/Rfam run as ncRNA annotation branches. Their outputs
are published under additional annotation outputs, recorded in provenance, used
by lncRNA candidate filtering, and summarized in ncRNA quality reporting.
They are not automatically merged into the final coding annotation.

Helixer, FLAIR, lncRNA/CPAT, SQANTI3, OMArk, eggNOG-mapper and InterProScan
are optional feature branches controlled by documented `run_*` parameters.
Optional branches should keep downstream reporting stable when disabled or
when no biological records are produced. Mikado and AEGIS are not optional
in the same sense: both default to enabled and the pipeline has no supported
path to a real final annotation with either disabled (see
[Mandatory Core Evidence](#mandatory-core-evidence)).

Mikado is the annotation consolidation step, not an alternative to AEGIS: it
receives the main annotation and transcript evidence families (Liftoff,
EGAPx, BRAKER3, STAR/StringTie, STAR/PsiCLASS, optional long-read StringTie,
optional FLAIR, optional Helixer) plus optional HISAT2/StringTie evidence
when `--run_hisat2 true`, then runs Mikado prepare, TransDecoder ORF
prediction, Mikado serialise and Mikado pick. Each source carries a numeric
`source_score` (and a `reference` flag for Liftoff/EGAPx) set in
`mikado_prepare`, which is what actually resolves competing gene models at a
locus — this replaced an earlier design where AEGIS merged the same sources
itself via whole-gene overlap-percentage exclusion, which had no way to keep
part of a gene model that didn't overlap. AEGIS now runs after Mikado, purely
to rename (`aegis_merge`, systematic Vitvi IDs) and then reconcile IDs against
the Liftoff-transferred previous annotation (`aegis_liftoff_gene_ids`, keeping
old IDs for genes with a confident, reciprocal-best-match correspondence via
`aegis overlap`). There is no AEGIS-vs-Mikado comparison anymore: AEGIS's
output is a renamed derivative of Mikado's, not an independent annotation.

## Staging And Publication

Processes should consume declared `path` inputs staged by Nextflow. Internal
directory scans of published output locations are avoided because they weaken
cache correctness and allow stale files to influence new runs.

`publishDir` is used to expose public results and maintain backward-compatible
filenames. It should not be used as a dependency mechanism between processes.

The main evidence path stages inputs for StringTie merges, GFFCompare,
BRAKER3, AGAT, Salmon and STAR/Minimap2 alignment stages. HISAT2 alignment
staging is present only when the optional `--run_hisat2 true` branch is enabled.

## Containers And Profiles

Runtime containers are centralized in `nextflow.config` and pinned with image
digests. Active modules should avoid Docker-specific mounts or options when
Nextflow staging can provide the required inputs. Profiles must remain usable
for local tests, Apptainer execution and Slurm production.

Resource policy belongs in config labels. Modules may express process labels,
but CPU, memory and cluster-specific tuning should be centralized whenever
possible.

## Validation And Test Data

Input validation runs before heavy computation. It checks required parameters,
input file existence, FASTA/GFF3 structure, samplesheet schemas, local FASTQ
and protein paths, minimal EGAPx YAML shape and constrained option values.

The `test` profile uses synthetic fixtures under `test-data/minimal` and runs
with the local executor in stub mode. It validates workflow structure, channel
shapes and file contracts. It does not validate biological correctness of
heavy external tools.

## Provenance

TITAN publishes provenance for the evidence bundle and workflow-level
versions. Tool-specific version capture should continue to expand as modules
are stabilized.

Public output names should stay deterministic and documented. When legacy
filenames are needed, prefer declared process outputs and `publishDir saveAs`
over manual copying into output directories.
