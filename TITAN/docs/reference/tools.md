# TITAN Tool Reference

This reference keeps tool-specific behavior out of the README while preserving
the operational details needed for production runs.

## Liftoff

Liftoff transfers the previous annotation from `--previous_assembly` and
`--previous_annotations` onto `--new_assembly`. TITAN publishes
`liftoff_previous_annotations.gff3` and `unmapped_features.txt`, then uses the
GFF3 as evidence for AEGIS, Mikado and Liftoff-derived splice-junction support
for FLAIR.

## EDTA

EDTA is mandatory. It annotates repetitive elements and produces the
hard-masked target genome consumed by AEGIS and downstream validation. TITAN
publishes the public masked genome as `assembly_masked.EDTA.fasta`.

## EGAPx

EGAPx is mandatory and is launched as a nested Nextflow run from the TITAN
`egapx` process. `--egapx_paramfile` points to an EGAPx YAML file; TITAN
validates at least `genome`, `taxid` and `organism`. EGAPx performs its own
masking and receives the original target assembly from the YAML, not the EDTA
masked assembly.

The nested runner writes `egapx_work` and `egapx_out` inside the task work
directory. Published outputs include genomic GFF3/GTF, proteins, CDS,
transcripts, ASN, the full `egapx_out/` directory and `versions.yml`.

## RNA-Seq Evidence

Short reads are prepared from local FASTQ files or ENA-resolved SRA accessions,
trimmed with fastp, aligned with STAR and HISAT2, and assembled with StringTie
and PsiCLASS. Salmon infers strandedness from Liftoff-derived CDS. Long-read
rows are detected from `library_layout=long`, aligned with Minimap2 and
assembled with StringTie.

Merged STAR/StringTie, HISAT2/StringTie, STAR/PsiCLASS and optional
Minimap2/StringTie GTFs are passed to AEGIS and Mikado. Per-sample alignments
and transcriptomes are published only when `--publish_intermediates true`.

## BRAKER3

BRAKER3 uses the masked genome, protein FASTA evidence and RNA-seq BAMs to
produce AUGUSTUS and GeneMark predictions. When long reads are present, TITAN
uses the long-read-aware BRAKER3 branch. Published outputs include
`augustus.hints.gff3`, `genemark.gtf`, `genemark_supported.gtf` and
`braker.gff3`.

## AEGIS

AEGIS is the primary final annotation integration step. It consumes named
evidence channels from Liftoff, EGAPx, BRAKER3, transcript assemblies, optional
long-read evidence, optional Helixer and optional FLAIR. It writes
`final_annotation.gff3`, `final_annotation_proteins_all.fasta` and
`final_annotation_proteins_main.fasta` under `${output_dir}/aegis_outputs`.

## Diamond2GO

Diamond2GO runs on the final AEGIS protein FASTAs and is part of the default
functional annotation path. Its outputs are published under
`${output_dir}/Diamond2GO_outputs`.

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

## Helixer

Helixer is optional (`--run_helixer true`) and predicts genes directly from the
EDTA soft-masked genome. Its GFF3 is published under
`${output_dir}/additional_annotations/helixer` and passed to AEGIS and Mikado
as optional evidence.

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

## tRNAscan-SE

tRNAscan-SE is optional (`--run_trnascan true`) and runs directly on the target
genome in parallel with evidence generation. TITAN runs eukaryotic mode, keeps
raw table/structure/isotype/statistics files, and converts the raw table to
`trna.gff3` with `scripts/trnascan_to_gff3.py`.

Outputs are published under `${output_dir}/additional_annotations/ncrna/trna`.
The tRNA GFF3 is used by lncRNA filtering and ncRNA QC, but is not
automatically merged into the AEGIS coding annotation.

## Infernal/Rfam ncRNA

Infernal/Rfam is optional (`--run_rfam true`) and runs directly on the target
genome. Stage Rfam offline once (`Rfam.cm`, `Rfam.clanin` and `cmpress`
indexes) and pass `--rfam_data_dir /absolute/path/to/rfam_data`.

TITAN splits the target FASTA by sequence, runs `cmsearch --cut_ga --rfam
--nohmmonly` independently on each split, merges all `rfam_hits.tbl` fragments,
and converts them once to `rfam_ncrna.gff3` with
`scripts/rfam_tblout_to_gff3.py`. Rfam outputs are published under
`${output_dir}/additional_annotations/ncrna/rfam`, used by lncRNA filtering and
ncRNA QC, and are not automatically merged into AEGIS.

## lncRNA Candidates

The lncRNA branch is optional (`--run_lncrna true`). It builds preliminary
candidate lncRNAs from merged transcript evidence after AEGIS, tRNAscan-SE and
Infernal/Rfam complete. Candidates are filtered by minimum length, excluded
when they overlap coding CDS, tRNA or Rfam intervals, and scored with the
bundled Plant-LncPipe CPAT-plant model.

Outputs include `lncrna_candidates.gff3`, `lncrna_candidates.gtf`,
`lncrna_candidates.fasta`, CPAT TSV/log files and summary TSV files under
`${output_dir}/additional_annotations/ncrna/lncrna`. This is a candidate
layer, not a final lncRNA annotation.

## Mikado Final Annotation Source

Mikado is optional (`--run_mikado true`) and produces an alternative final GFF3
source in parallel with AEGIS. It receives the same major evidence families as
AEGIS: Liftoff, EGAPx, BRAKER3, STAR/StringTie, HISAT2/StringTie,
STAR/PsiCLASS, optional long-read StringTie, optional FLAIR and optional
Helixer.

The graph runs Mikado configure/prepare, TransDecoder LongOrfs/Predict when
`--run_transdecoder true`, then Mikado serialise and pick. Outputs are
published under `${output_dir}/final_annotations/mikado`. AEGIS and Mikado are
separate final GFF3 sources; TITAN adds an AEGIS-vs-Mikado overlap summary to
MultiQC when Mikado is enabled.

## FLAIR Long-Read Isoforms

FLAIR is optional (`--run_flair true`) and runs only when long-read samples are
present. It uses Liftoff as splice-junction correction evidence, publishes
per-sample and merged isoform GTF/FASTA files under
`${output_dir}/additional_annotations/flair`, and passes merged isoforms to
AEGIS and Mikado as optional transcript evidence.

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

Outputs are published under `${output_dir}/additional_annotations/sqanti3` and
`${output_dir}/quality_report/sqanti3`.

## Quality Report

TITAN closes with quality reporting over the final AEGIS annotation and enabled
optional branches:

* BUSCO protein-mode completeness when `--run_busco true`;
* OMArk consistency/contamination checks when `--run_omark true`;
* AGAT structural statistics on `final_annotation.gff3`;
* ncRNA counts from tRNAscan-SE and Rfam outputs;
* expression support validation with Salmon when `--run_expression_validation true`;
* SQANTI3 summaries when enabled and long-read evidence exists;
* AEGIS-vs-Mikado comparison when Mikado is enabled;
* final MultiQC HTML aggregation.

BUSCO needs an offline lineage dataset and currently has no
`--prepare-busco-data` launcher flag. OMArk requires an offline OMAmer database
at `--omark_data_dir` containing `omamer.h5`; the launchers can prepare it with
`--prepare-omark-data`.
