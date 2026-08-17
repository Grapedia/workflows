# TITAN input and output tree layouts

## Table Of Contents

- [Input Tree](#input-tree)
- [Output Tree](#output-tree)
- [Primary Annotation](#primary-annotation)
- [High-Confidence Monoexonic-Filtered Annotation](#high-confidence-monoexonic-filtered-annotation)
- [Quality, Validation And Provenance](#quality-validation-and-provenance)
- [Functional Annotation](#functional-annotation)
- [Additional Annotations](#additional-annotations)
- [Evidence: Core Tracks](#evidence-core-tracks)
- [Evidence: EGAPx](#evidence-egapx)
- [Evidence: Mikado Consolidated Annotation](#evidence-mikado-consolidated-annotation)
- [Run Info: Provenance And Runtime Outputs](#run-info-provenance-and-runtime-outputs)


This page gives a file-tree view of the inputs TITAN expects and the outputs it
publishes. Paths are examples; use absolute paths for production and HPC runs.

## Input Tree

Recommended project layout:

```text
project/
  assemblies/
    target.fa
    previous.fa
  annotations/
    previous.gff3
  rnaseq/
    RNAseq_samplesheet.csv
    leaf_single.fastq.gz
    berry_paired_1.fastq.gz
    berry_paired_2.fastq.gz
    isoseq_leaf.fastq.gz
  proteins/
    protein_samplesheet.csv
    araport.fa
    swissprot_plants.fa
  egapx/
    input_egapx.yaml
  titan_out/
```

`assemblies/target.fa` is `--new_assembly`, the genome to annotate.
`assemblies/previous.fa` is `--previous_assembly`, the genome associated with
the transferred annotation. `annotations/previous.gff3` is
`--previous_annotations` and must match the previous assembly sequence IDs.

`rnaseq/RNAseq_samplesheet.csv` is `--RNAseq_samplesheet`. Local files are
found under `--RNAseq_data_dir` from each `sampleID`: `<sampleID>.fastq.gz` for
single-end short reads, `<sampleID>_1.fastq.gz` and `<sampleID>_2.fastq.gz` for
paired-end short reads, and `<sampleID>.fastq.gz` or `<sampleID>.fasta` for
long-read transcript evidence.

`proteins/protein_samplesheet.csv` is `--protein_samplesheet`; each `filename`
entry points to a protein FASTA used by BRAKER3/ProtHint. `egapx/input_egapx.yaml`
is `--egapx_paramfile` and must contain at least `genome`, `taxid` and
`organism`.

Optional reference caches can be kept outside `output_dir`, for example:

```text
project/
  .egapx_runner/
  .egapx_cache/
  .eggnog_data/
  .helixer_models/
  .interproscan_data/
  .rfam_data/
  .omark_data/
```

These directories are selected with launcher preparation flags or native
Nextflow parameters. They are inputs/reference data, not TITAN results.

## Output Tree

The exact tree depends on enabled optional branches and
`--publish_intermediates`. A complete production run publishes the main result
families below under `--output_dir`, in 5 numbered top-level folders — sorted
in reading order, so `ls` lists them in the order you should actually look at
them.

```text
${output_dir}/
  01_final_annotation/
    primary/
      final_annotation.gff3
      final_annotation_proteins_all.fasta
      final_annotation_proteins_main.fasta
      liftoff_gene_id_correspondence.tsv
      versions.yml
    high_confidence_monoexonic/
      final_annotation.high_confidence.gff3
      final_annotation_proteins_main.high_confidence.fasta
      monoexonic_gene_confidence.tsv
      monoexonic_gene_confidence_summary.json
      versions.yml
      agat_stats/
      busco/
    quality_report/
      titan_multiqc_report.html
      titan_multiqc_report_data/
      agat_stats/
      busco/
      expression_validation/
      ncrna_annotations/
      omark/
      sqanti3/
      versions.yml
    validation/
      final_annotation_validation.json
      final_annotation_validation.txt
      versions.yml
  02_functional_annotation/
    diamond2go/
    eggnog/
    interproscan/
  03_additional_annotations/
    ncrna/
    flair/
    helixer/
    sqanti3/
  04_evidence/
    assembly_masked.EDTA.fasta
    liftoff_previous_annotations.gff3
    unmapped_features.txt
    versions.yml
    gene_prediction/
      augustus.hints.gff3
      genemark.gtf
      genemark_supported.gtf
      braker.gff3
      versions.yml
    transcript_assemblies/
      merged_star_stringtie_stranded_default.gtf
      merged_star_stringtie_stranded_alt.gtf
      merged_star_stringtie_unstranded_default.gtf
      merged_star_stringtie_unstranded_alt.gtf
      merged_star_psiclass_stranded.gtf
      merged_star_psiclass_unstranded.gtf
      merged_minimap2_stringtie_long_reads_default.gtf
      merged_minimap2_stringtie_long_reads_alt.gtf
      # when --run_hisat2 true:
      merged_hisat2_stringtie_stranded_default.gtf
      merged_hisat2_stringtie_stranded_alt.gtf
      merged_hisat2_stringtie_unstranded_default.gtf
      merged_hisat2_stringtie_unstranded_alt.gtf
      versions_star_psiclass.yml
      versions_star_stringtie.yml
      versions_hisat2_stringtie.yml
      versions_minimap2_stringtie_long_reads.yml
    egapx/
    mikado/
  05_run_info/
    provenance/
    intermediate_files/
    nextflow_reports/
  versions.yml
```

`01_final_annotation/primary/` and `01_final_annotation/high_confidence_monoexonic/`
are the two final annotation candidates to compare (see
[High-Confidence Monoexonic-Filtered Annotation](#high-confidence-monoexonic-filtered-annotation)),
with their shared `quality_report/` and `validation/` right alongside them —
so `01_final_annotation/` as a whole is what to check first on a completed
run. `02_functional_annotation/` and `03_additional_annotations/` come next.
`04_evidence/` groups everything that feeds Mikado/AEGIS but is not itself a
final annotation — including the EGAPx run and the pre-AEGIS Mikado
consolidation (see [Evidence: EGAPx](#evidence-egapx) and
[Evidence: Mikado Consolidated Annotation](#evidence-mikado-consolidated-annotation)).
`05_run_info/` is provenance/debug material, not required reading. Some
unstranded or long-read merged GTF files can be empty when the corresponding
sample class is absent.

## Primary Annotation

```text
${output_dir}/01_final_annotation/primary/
  final_annotation.gff3
  final_annotation_proteins_all.fasta
  final_annotation_proteins_main.fasta
  liftoff_gene_id_correspondence.tsv
  versions.yml
```

This is the primary output family. `final_annotation.gff3` is Mikado's
consolidated gene set, renamed and tidied by AEGIS (systematic `Vitvi` IDs,
with old gene IDs carried over from the previous annotation via Liftoff
where a confident correspondence exists). The protein FASTAs contain all
translated proteins and the main representative protein set.
`liftoff_gene_id_correspondence.tsv` records, per gene, whether its ID was
carried over (`carried_over_from_liftoff`, with the match's overlap score)
or freshly assigned (`new_vitvi_id`).

## High-Confidence Monoexonic-Filtered Annotation

```text
${output_dir}/01_final_annotation/high_confidence_monoexonic/
  final_annotation.high_confidence.gff3
  final_annotation_proteins_main.high_confidence.fasta
  monoexonic_gene_confidence.tsv
  monoexonic_gene_confidence_summary.json
  versions.yml
  agat_stats/
    agat_stats.txt
    versions.yml
  busco/
    busco_short_summary.txt
    busco_full_table.tsv
    busco_missing_busco_list.tsv
    versions.yml
```

This is the second final-annotation candidate:
`01_final_annotation/primary/final_annotation.gff3` with every unsupported
single-exon gene removed. `monoexonic_gene_confidence.tsv` classifies every
single-exon gene in the primary annotation by supporting evidence
(expression, functional domains, Liftoff conservation, BUSCO orthology, TE
overlap) without touching the annotation itself; see
[tools reference](../reference/tools.md#single-exon-gene-confidence).
`monoexonic_gene_confidence_summary.json` is the run-level summary.

`agat_stats/` and `busco/` here run on this filtered GFF3/protein set, so they
can be compared directly against
`01_final_annotation/quality_report/agat_stats/` and
`01_final_annotation/quality_report/busco/` (the primary candidate's) to
decide which annotation to keep or publish.

## Quality, Validation And Provenance

```text
${output_dir}/01_final_annotation/quality_report/
  titan_multiqc_report.html
  titan_multiqc_report_data/
  agat_stats/
  busco/
  expression_validation/
  ncrna_annotations/
  omark/
  sqanti3/
  versions.yml
${output_dir}/01_final_annotation/validation/
  final_annotation_validation.json
  final_annotation_validation.txt
  versions.yml
```

`quality_report/titan_multiqc_report.html` is the main QC entry point;
`quality_report/agat_stats/` and `quality_report/busco/` assess the primary
`01_final_annotation/primary/` candidate specifically. See
[High-Confidence Monoexonic-Filtered Annotation](#high-confidence-monoexonic-filtered-annotation)
for the equivalent AGAT/BUSCO results on the second candidate. `validation/`
reports structural checks on the final annotation. Input/evidence provenance
manifests live under `05_run_info/provenance/` (see
[Run Info: Provenance And Runtime Outputs](#run-info-provenance-and-runtime-outputs)).

## Functional Annotation

```text
${output_dir}/02_functional_annotation/diamond2go/
  final_annotation_proteins_all.diamond2go.tsv
  final_annotation_proteins_main.diamond2go.tsv
  versions.yml
${output_dir}/02_functional_annotation/eggnog/
  final_annotation_proteins_all.emapper.annotations
  final_annotation_proteins_main.emapper.annotations
  final_annotation_proteins_all.emapper.seed_orthologs
  final_annotation_proteins_main.emapper.seed_orthologs
  final_annotation_proteins_all.emapper.orthologs
  final_annotation_proteins_main.emapper.orthologs
  final_annotation_proteins_all.emapper.annotations.xlsx
  final_annotation_proteins_main.emapper.annotations.xlsx
  versions.yml
${output_dir}/02_functional_annotation/interproscan/
  final_annotation_proteins_all.tsv
  final_annotation_proteins_main.tsv
  final_annotation_proteins_all.gff3
  final_annotation_proteins_main.gff3
  final_annotation_proteins_all.json
  final_annotation_proteins_main.json
  versions.yml
```

All three annotate `01_final_annotation/primary/`'s proteins (the primary
candidate, not the high-confidence-filtered one). Each tool runs on the
protein FASTA, so its raw output keys every row on the FASTA header
(`<gene_id>_t<NNN>_CDS<N>.prot`, a CDS-record ID, e.g.
`Vitvichr00g00010_t001_CDS1.prot`); every one of these files is rewritten at
the end of its job to collapse that ID back to the bare gene ID
(`Vitvichr00g00010`) everywhere it appears — TSV query columns, the
InterProScan GFF3's seqid/`ID=`/`##sequence-region`, its JSON `xref`
entries, and the eggNOG-mapper `.xlsx` exports — so results join directly
onto `final_annotation.gff3` by gene ID. `_all` files can have several rows
per gene when the primary annotation has multiple isoforms for it.
Diamond2GO is part of the default graph. eggNOG-mapper and InterProScan are
optional and only appear when enabled.

## Additional Annotations

```text
${output_dir}/03_additional_annotations/
  ncrna/
    trna/
      trna.gff3
      trnascan.out
      trnascan.stats
      trnascan.struct
      trnascan.isotype
      versions.yml
    rfam/
      rfam_hits.tbl
      rfam_search.out
      rfam_ncrna.gff3
      versions.yml
    lncrna/
      lncrna_candidates.gff3
      lncrna_candidates.gtf
      lncrna_candidates.fasta
      lncrna_classification_summary.tsv
      lncrna_candidates_mqc.tsv
      cpat_plant.output.ORF_prob.tsv
      cpat_plant.output.ORF_prob.best.tsv
      cpat_plant.output.no_ORF.txt
      CPAT_run_info.log
      versions.yml
  flair/
    flair_isoforms.gtf
    flair_isoforms.fa
    samples/
    versions.yml
  helixer/
    helixer.gff3
    versions.yml
  sqanti3/
    <source_label>/
```

These branches are optional except where enabled by the run configuration.
Rfam and tRNAscan-SE feed ncRNA summaries and lncRNA filtering. FLAIR and
SQANTI3 require long-read evidence and their own optional settings. Helixer and
lncRNA outputs are only present when the corresponding branches are enabled.

## Evidence: Core Tracks

```text
${output_dir}/04_evidence/
  assembly_masked.EDTA.fasta
  liftoff_previous_annotations.gff3
  unmapped_features.txt
  gene_prediction/
    augustus.hints.gff3
    genemark.gtf
    genemark_supported.gtf
    braker.gff3
```

`assembly_masked.EDTA.fasta` is the EDTA-masked assembly. Liftoff outputs are
the transferred previous annotation and the list of unmapped features. BRAKER3
publishes AUGUSTUS and GeneMark evidence under `04_evidence/gene_prediction/`
for direct inspection.

Merged transcript evidence is published under `04_evidence/transcript_assemblies/`:

```text
${output_dir}/04_evidence/transcript_assemblies/
  merged_star_stringtie_*.gtf
  merged_star_psiclass_*.gtf
  merged_minimap2_stringtie_long_reads_*.gtf
  # when --run_hisat2 true:
  merged_hisat2_stringtie_*.gtf
```

These files summarize short-read STAR/PsiCLASS/StringTie evidence and long-read
Minimap2/StringTie evidence. HISAT2/StringTie merged tracks are published only
when the optional branch is enabled with `--run_hisat2 true`. This directory
is the single copy of these merges; earlier TITAN versions also published a
byte-identical, pre-rename copy under `${output_dir}/tmp/`, which has been
removed.

## Evidence: EGAPx

```text
${output_dir}/04_evidence/egapx/
  egapx.complete.genomic.gff3
  egapx.complete.genomic.gtf
  egapx.complete.proteins.faa
  egapx.complete.cds.fna
  egapx.complete.transcripts.fna
  egapx.annotated_genome.asn
  egapx_out/
  versions.yml
```

The `egapx.complete.*` files are copied from the nested EGAPx run into stable
TITAN filenames. `egapx_out/` is the nested EGAPx output directory. This is a
standalone NCBI EGAPx annotation run, kept as evidence/comparison material
alongside the Mikado/AEGIS evidence tracks, not a TITAN final annotation.

## Evidence: Mikado Consolidated Annotation

```text
${output_dir}/04_evidence/mikado/
  final_mikado_annotation.gff3
  mikado.loci.gff3
  mikado.subloci.gff3
  versions.yml
  intermediate/
  transdecoder/
```

Mikado consolidates every evidence source into one non-redundant, per-locus-
scored gene set; `final_mikado_annotation.gff3` is what AEGIS then renames
(Vitvi IDs, with liftoff ID carryover) and tidies into
`01_final_annotation/primary/final_annotation.gff3`. This branch defaults to
enabled (`--run_mikado true`) and is effectively required — disabling it
leaves `aegis_merge` nothing to rename and the run fails there. It is
evidence for `01_final_annotation/primary/`, not a final annotation in its
own right.

## Run Info: Provenance And Runtime Outputs

```text
${output_dir}/05_run_info/
  provenance/
    evidence_manifest.json
    additional_annotations_manifest.json
    versions.yml
  intermediate_files/
    aegis/                  # aegis_merge's Vitvi-only pass (vitvi_only_*) and
                             # the aegis overlap correspondence table
    braker3/
    evidence_data/
      EDTA/
      RNAseq_alignments/
        STAR/
        minimap2/
        HISAT2/              # when --run_hisat2 true
      RNAseq_data/
        trimmed_data/
      transcriptomes/
        STAR_PsiCLASS/
        StringTie/
      star_databases/
      hisat2_databases/      # when --run_hisat2 true
    liftoff/
    salmon_strand/
  nextflow_reports/
    <run_name>.dag.html
    progress.log
```

Provenance manifests record inputs, selected outputs, versions and checksums.
`intermediate_files/` is controlled by `--publish_intermediates`. It is
useful for debugging, manual inspection and downstream reuse, but the main
deliverables are `01_final_annotation/`, `02_functional_annotation/` and
`05_run_info/provenance/` — `04_evidence/` (see
[Evidence: Core Tracks](#evidence-core-tracks)) and
`05_run_info/intermediate_files/` sit one step below those as supporting
material worth keeping but not required reading.

The Nextflow `work/` directory and `.nextflow.log` are runtime state, not
published biological outputs. Keep them for resume/debugging, but do not treat
them as deliverables.
