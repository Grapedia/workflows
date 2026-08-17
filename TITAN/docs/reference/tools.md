# TITAN Tool Reference

## Table Of Contents

- [Input Validation](#input-validation)
- [tRNAscan-SE](#trnascan-se)
- [Infernal/Rfam ncRNA](#infernalrfam-ncrna)
- [RNA-Seq Evidence](#rna-seq-evidence)
- [Liftoff](#liftoff)
- [EGAPx](#egapx)
- [FLAIR Long-Read Isoforms](#flair-long-read-isoforms)
- [EDTA](#edta)
- [BRAKER3](#braker3)
- [Helixer](#helixer)
- [Mikado Final Annotation Source](#mikado-final-annotation-source)
- [AEGIS](#aegis)
- [lncRNA Candidates](#lncrna-candidates)
- [SQANTI3 Long-Read Isoform QC](#sqanti3-long-read-isoform-qc)
- [Diamond2GO](#diamond2go)
- [eggNOG-mapper](#eggnog-mapper)
- [InterProScan](#interproscan)
- [Quality Report](#quality-report)


This reference keeps tool-specific behavior out of the README while preserving
the operational details needed for production runs. Sections follow the
pipeline orchestration order in `workflows/titan.nf`: validation, early
genome-wide optional branches, evidence generation, final integration,
post-merge annotations, functional annotation and quality reporting.

For the exact command lines launched by each process, see
[tool-commands.md](tool-commands.md).

## Input Validation

TITAN starts by rejecting deprecated public modes, checking required parameter
presence, validating input file paths, parsing the RNA-seq and protein
samplesheets, and checking optional external-data directories before heavy
work begins. This includes fail-fast checks for eggNOG-mapper, Helixer,
InterProScan and CPAT/lncRNA inputs when their branches are enabled.

## tRNAscan-SE

tRNAscan-SE is optional (`--run_trnascan true`) and runs directly on the target
genome near the start of the graph. TITAN runs eukaryotic mode, keeps raw
table/structure/isotype/statistics files, and converts the raw table to
`trna.gff3` with `scripts/trnascan_to_gff3.py`.

Outputs are published under `${output_dir}/03_additional_annotations/ncrna/trna`.
The tRNA GFF3 is used by lncRNA filtering and ncRNA QC, but is not
automatically merged into the final coding annotation.

## Infernal/Rfam ncRNA

Infernal/Rfam is optional (`--run_rfam true`) and runs directly on the target
genome near the start of the graph. Stage Rfam offline once (`Rfam.cm`,
`Rfam.clanin` and `cmpress` indexes) and pass
`--rfam_data_dir /absolute/path/to/rfam_data`.

TITAN splits the target FASTA by sequence, runs `cmsearch --cut_ga --rfam
--nohmmonly` independently on each split, merges all `rfam_hits.tbl` fragments,
and converts them once to `rfam_ncrna.gff3` with
`scripts/rfam_tblout_to_gff3.py`. Rfam outputs are published under
`${output_dir}/03_additional_annotations/ncrna/rfam`, used by lncRNA filtering and
ncRNA QC, and are not automatically merged into the final coding annotation.

## RNA-Seq Evidence

Short reads are prepared from local FASTQ files or ENA-resolved SRA accessions,
trimmed with fastp, aligned with STAR, and assembled with StringTie and
PsiCLASS. Salmon infers strandedness from Liftoff-derived CDS. HISAT2/StringTie
is an optional short-read evidence branch disabled by default; enable it with
`--run_hisat2 true`. Long-read rows are detected from `library_layout=long`,
aligned with Minimap2 and assembled with StringTie.

Merged STAR/StringTie, STAR/PsiCLASS and optional Minimap2/StringTie GTFs are
passed to Mikado (AEGIS no longer consumes raw evidence directly — see
[AEGIS](#aegis)). When enabled, HISAT2/StringTie also feeds Mikado, lncRNA
candidates and provenance. Per-sample alignments and transcriptomes are
published only when `--publish_intermediates true`.

## Liftoff

Liftoff transfers the previous annotation from `--previous_assembly` and
`--previous_annotations` onto `--new_assembly`. TITAN publishes
`liftoff_previous_annotations.gff3` and `unmapped_features.txt`, then uses the
GFF3 as Mikado consolidation evidence, as the correspondence source for
AEGIS's liftoff gene-ID carryover (see [AEGIS](#aegis)), and as
Liftoff-derived splice-junction support for FLAIR.

## EGAPx

EGAPx is mandatory and is launched as a nested Nextflow run from the TITAN
`egapx` process. `--egapx_paramfile` points to an EGAPx YAML file; TITAN
validates at least `genome`, `taxid` and `organism`. EGAPx performs its own
masking and receives the original target assembly from the YAML, not the EDTA
masked assembly.

The nested runner writes `egapx_work` and `egapx_out` inside the task work
directory. Published outputs include genomic GFF3/GTF, proteins, CDS,
transcripts, ASN, the full `egapx_out/` directory and `versions.yml`.

## FLAIR Long-Read Isoforms

FLAIR is optional (`--run_flair true`) and runs only when long-read samples are
present. It uses Liftoff as splice-junction correction evidence, publishes
per-sample and merged isoform GTF/FASTA files under
`${output_dir}/03_additional_annotations/flair`, and passes merged isoforms to
Mikado as optional transcript evidence.

## EDTA

EDTA is mandatory. It annotates repetitive elements and produces the
hard-masked target genome consumed by Mikado, AEGIS and downstream
validation. TITAN publishes the public masked genome as
`assembly_masked.EDTA.fasta`.

## BRAKER3

BRAKER3 uses the target genome, protein FASTA evidence and RNA-seq BAMs to
produce AUGUSTUS and GeneMark predictions. When long reads are present, TITAN
uses the long-read-aware BRAKER3 branch. Published outputs include
`augustus.hints.gff3`, `genemark.gtf`, `genemark_supported.gtf` and
`braker.gff3`.

## Helixer

Helixer is optional (`--run_helixer true`) and predicts genes directly from the
EDTA soft-masked genome. Its GFF3 is published under
`${output_dir}/03_additional_annotations/helixer` and passed to Mikado as
optional evidence (Mikado is the only step that consumes raw evidence
sources; AEGIS only sees Mikado's consolidated output — see
[AEGIS](#aegis) and [Mikado](#mikado-final-annotation-source)).

Use `--helixer_model_dir` with a staged lineage model and optionally
`--helixer_model` (`vertebrate`, `land_plant`, `fungi` or `invertebrate`;
default `land_plant`). CPU is the default; `--helixer_use_gpu true` requests
GPU runtime options and requires a GPU-visible node.

Fetch a model once with:

```bash
scripts/download_helixer_model.sh \
  --model-dir /absolute/path/to/helixer_models \
  --container <container_helixer> \
  --lineage land_plant
```

## Mikado Final Annotation Source

Mikado is the annotation consolidation step: it reconciles every evidence
source into one non-redundant, per-locus-scored gene set, which AEGIS then
renames and tidies (see [AEGIS](#aegis) below). It defaults to enabled
(`--run_mikado true`) and is effectively required — if disabled, `mikado_pick`
produces an empty annotation and `aegis_merge` fails fast rather than silently
publishing nothing. It receives Liftoff, EGAPx, BRAKER3, STAR/StringTie,
STAR/PsiCLASS, optional long-read StringTie, optional FLAIR and optional
Helixer. If `--run_hisat2 true`, TITAN also adds HISAT2/StringTie sources to
Mikado; otherwise skipped HISAT2 records are written to provenance.

Each source is given a calibrated numeric priority (`source_score` in
`modules/mikado.nf`'s `mikado_prepare`), highest first: EGAPx (20) > Liftoff
(19) > BRAKER-AUGUSTUS (18) > BRAKER-GeneMark (17) > Helixer (16) > stranded
assemblies (STAR/StringTie default 15, alt 14, STAR/PsiCLASS 13, then their
optional HISAT2/StringTie counterparts 12/11) > FLAIR (10) > raw long reads
(default 9, alt 8) > unstranded assemblies (STAR/PsiCLASS 7, STAR/StringTie
default 6, alt 5, then optional HISAT2/StringTie 4/3). Liftoff/EGAPx are
additionally flagged as `reference` quality (excluded only for outright
mistakes, not overlap competition). This scoring is what determines which
transcript wins at a given locus, replacing AEGIS's older whole-gene
overlap-threshold merge.

The graph runs Mikado configure/prepare, TransDecoder LongOrfs/Predict when
`--run_transdecoder true` (also effectively required — without ORFs, Mikado's
CDS-based scoring has nothing to score and rejects most transcripts), then
Mikado serialise and pick. Outputs are published as evidence under
`${output_dir}/04_evidence/mikado`; `final_mikado_annotation.gff3` is the file
AEGIS renames into the final annotation.

## AEGIS

AEGIS no longer merges evidence directly; Mikado does that (see
[Mikado Final Annotation Source](#mikado-final-annotation-source) above).
AEGIS runs in two steps on Mikado's consolidated output:

1. **`aegis_merge`** (rename/tidy) — runs `aegis rename --prefix <aegis_gene_id_prefix>`
   on `final_mikado_annotation.gff3` to assign systematic, position-based gene
   IDs (default prefix `Vitvi`, independent of Mikado's own IDs), then
   `aegis tidy --standard-features`. This intermediate, Vitvi-only annotation
   is not published by default (only under
   `${output_dir}/05_run_info/intermediate_files/aegis` with `--publish_intermediates true`,
   as `vitvi_only_final_annotation.gff3`).
2. **`aegis_liftoff_gene_ids`** (ID carryover) — runs `aegis overlap` between
   that Vitvi-renamed annotation and the Liftoff-transferred previous
   annotation (`evidence_data.liftoff_gff3`), with `previous_annotations` as
   `--original-annotation-files` so AEGIS can also check synteny conservation
   against the true old-assembly gene order. For every gene with a confident,
   reciprocal-best-match correspondence to an old gene (AEGIS's default
   `overlap_score` threshold), the old gene ID replaces the fresh Vitvi ID
   (cascaded to transcript/exon/CDS/UTR IDs and protein FASTA headers via
   `scripts/apply_liftoff_gene_ids.py`); genes without a confident match keep
   their Vitvi ID. Ambiguous matches (e.g. a gene fusion/split) are treated as
   no-match and also keep a fresh Vitvi ID, rather than risk a misleading old
   ID assignment.

This step writes the files that are actually published under
`${output_dir}/01_final_annotation/primary`: `final_annotation.gff3`,
`final_annotation_proteins_all.fasta`, `final_annotation_proteins_main.fasta`
and `liftoff_gene_id_correspondence.tsv` (per-gene decision and match score,
for provenance).

## lncRNA Candidates

The lncRNA branch is optional (`--run_lncrna true`). It builds preliminary
candidate lncRNAs from merged transcript evidence after AEGIS, tRNAscan-SE and
Infernal/Rfam complete. Candidates are filtered by minimum length, excluded
when they overlap coding CDS, tRNA or Rfam intervals, and scored with the
bundled Plant-LncPipe CPAT-plant model.

Outputs include `lncrna_candidates.gff3`, `lncrna_candidates.gtf`,
`lncrna_candidates.fasta`, CPAT TSV/log files and summary TSV files under
`${output_dir}/03_additional_annotations/ncrna/lncrna`. This is a candidate
layer, not a final lncRNA annotation.

## SQANTI3 Long-Read Isoform QC

SQANTI3 is optional (`--run_sqanti3 true`) and runs after AEGIS. It evaluates
StringTie/Minimap2 long-read transcripts and FLAIR isoforms against the final
AEGIS annotation and target genome. When long reads are absent, or when one
source has no isoforms, TITAN emits zero-count sentinel summaries so MultiQC
remains stable.

The pinned SQANTI3 image needs a task-local `libbz2.so.1` compatibility
symlink for `gtfToGenePred`; TITAN creates it from the container-visible path
configured by `--sqanti3_libbz2_path` (default
`/usr/local/lib/libbz2.so.1.0.8`). This is deliberately not a host bind path.
Set the parameter to `false` only when a replacement SQANTI3 image already
ships `libbz2.so.1`.

Outputs are published under `${output_dir}/03_additional_annotations/sqanti3` and
`${output_dir}/01_final_annotation/quality_report/sqanti3`.

## Diamond2GO

Diamond2GO runs on the final AEGIS protein FASTAs and is part of the default
functional annotation path. Its outputs are published under
`${output_dir}/02_functional_annotation/diamond2go`.

Diamond2GO, eggNOG-mapper and InterProScan (below) all run on the protein
FASTAs, so each tool's raw output keys every row on the FASTA header
(`<gene_id>_t<NNN>_CDS<N>.prot`, a CDS-record ID), not the gene ID. Each
module's script rewrites its own outputs in place at the end of the job
(`sed`, plus a Python step in `eggnog_mapper` to regenerate the `.xlsx`
exports from the now-cleaned annotations TSV) to collapse that ID back to
the bare gene ID everywhere it appears, so results join directly onto
`final_annotation.gff3` by gene ID.

## eggNOG-mapper

eggNOG-mapper is optional (`--run_eggnog_mapper true`) and runs on the AEGIS
protein FASTAs in parallel with Diamond2GO. It requires a pre-downloaded
database directory passed with `--eggnog_data_dir`. Sensitivity and optional
taxonomic scope are controlled by `--eggnog_mapper_sensmode` and
`--eggnog_mapper_tax_scope`.

Fetch the database once with:

```bash
scripts/download_eggnog_data.sh --data-dir /absolute/path/to/eggnog_data
```

The launchers can run the same preparation with `--prepare-eggnog-data`.

## InterProScan

InterProScan is optional (`--run_interproscan true`) and runs on the AEGIS
protein FASTAs. It requires pre-downloaded member database data passed with
`--interproscan_data_dir`. TITAN uses offline mode with GO terms and pathways
enabled and disables the online precalculated match lookup service.

Fetch the member database once with:

```bash
scripts/download_interproscan_data.sh --data-dir /absolute/path/to/interproscan_data
```

The launchers can run the same preparation with
`--prepare-interproscan-data`.

## Quality Report

TITAN closes with quality reporting over the final AEGIS annotation and enabled
optional branches:

* BUSCO protein-mode completeness when `--run_busco true`;
* OMArk consistency/contamination checks when `--run_omark true`;
* AGAT structural statistics on `final_annotation.gff3`;
* ncRNA counts from tRNAscan-SE and Rfam outputs;
* expression support validation with Salmon when `--run_expression_validation true`;
* SQANTI3 summaries when enabled and long-read evidence exists;
* final MultiQC HTML aggregation.

BUSCO needs an offline lineage dataset at `--busco_data_dir`; the launchers can
prepare it with `--prepare-busco-data`, which uses the pinned BUSCO container.
OMArk requires an offline OMAmer database at `--omark_data_dir` containing
`omamer.h5`; the launchers can prepare it with `--prepare-omark-data`.

## Single-exon gene confidence

Ab initio predictors (BRAKER/Helixer) are known to over-predict single-exon
genes in plant genomes; a large fraction of `final_annotation.gff3`'s genes
can be monoexonic without independent support. `flag_low_confidence_monoexonic_genes`
(`scripts/flag_low_confidence_monoexonic_genes.py`) cross-references every
single-exon gene against evidence TITAN already computes elsewhere in the
pipeline:

* Salmon expression support (`expression_support_summary`'s unsupported gene list);
* functional domain/GO hits (InterProScan, eggNOG-mapper, Diamond2GO);
* gene IDs conserved from the previous annotation via liftoff;
* BUSCO single-copy/duplicated/fragmented orthology;
* overlap with an EDTA-annotated transposable element (negative signal only).

Each single-exon gene is classified into one of three tiers and reported in
`${output_dir}/01_final_annotation/high_confidence_monoexonic/`:

* `supported` - at least one positive evidence signal; kept as-is;
* `unsupported_te_overlap` - no positive evidence, overlaps a TE;
* `unsupported_grey_zone` - no positive evidence, no TE overlap either.

`final_annotation.gff3` itself is never modified. The process also writes a
"high confidence" variant with every non-`supported` single-exon gene (and
its child features) removed:
`final_annotation.high_confidence.gff3` /
`final_annotation_proteins_main.high_confidence.fasta`. BUSCO and AGAT stats
are run on this filtered variant too
(`busco_high_confidence_monoexonic`, `agat_stats_high_confidence_monoexonic`,
published under `${output_dir}/01_final_annotation/high_confidence_monoexonic/busco`
and `.../01_final_annotation/high_confidence_monoexonic/agat_stats`), so
completeness and gene-set structure can be compared side by side against the
unfiltered annotation's own `01_final_annotation/quality_report/busco`/
`01_final_annotation/quality_report/agat_stats`
results before deciding whether the filtered variant should replace the
primary annotation.
