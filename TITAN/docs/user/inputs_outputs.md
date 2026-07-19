# TITAN input and output tree layouts

## Table Of Contents

- [Input Tree](#input-tree)
- [Output Tree](#output-tree)
- [Primary Annotation](#primary-annotation)
- [Core Evidence Outputs](#core-evidence-outputs)
- [EGAPx Outputs](#egapx-outputs)
- [Additional Annotations](#additional-annotations)
- [Optional Final Annotation Source](#optional-final-annotation-source)
- [Functional Annotation](#functional-annotation)
- [Quality, Validation And Provenance](#quality-validation-and-provenance)
- [Intermediate And Runtime Outputs](#intermediate-and-runtime-outputs)


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
families below under `--output_dir`.

```text
${output_dir}/
  aegis_outputs/
    final_annotation.gff3
    final_annotation_proteins_all.fasta
    final_annotation_proteins_main.fasta
    versions.yml
  assembly_masked.EDTA.fasta
  liftoff_previous_annotations.gff3
  unmapped_features.txt
  augustus.hints.gff3
  genemark.gtf
  genemark_supported.gtf
  braker.gff3
  merged_hisat2_stringtie_stranded_default.gtf
  merged_hisat2_stringtie_stranded_alt.gtf
  merged_hisat2_stringtie_unstranded_default.gtf
  merged_hisat2_stringtie_unstranded_alt.gtf
  merged_star_stringtie_stranded_default.gtf
  merged_star_stringtie_stranded_alt.gtf
  merged_star_stringtie_unstranded_default.gtf
  merged_star_stringtie_unstranded_alt.gtf
  merged_star_psiclass_stranded.gtf
  merged_star_psiclass_unstranded.gtf
  merged_minimap2_stringtie_long_reads_default.gtf
  merged_minimap2_stringtie_long_reads_alt.gtf
  egapx/
  additional_annotations/
  final_annotations/
  Diamond2GO_outputs/
  EggNOG_outputs/
  InterProScan_outputs/
  quality_report/
  validation/
  provenance/
  intermediate_files/
  tmp/
  nextflow_reports/
  versions.yml
```

Top-level evidence tracks are intentionally kept near `aegis_outputs/` because
they are direct AEGIS inputs or easy-to-inspect evidence products. Some
unstranded or long-read merged GTF files can be empty when the corresponding
sample class is absent.

## Primary Annotation

```text
${output_dir}/aegis_outputs/
  final_annotation.gff3
  final_annotation_proteins_all.fasta
  final_annotation_proteins_main.fasta
  versions.yml
```

This is the primary output family. `final_annotation.gff3` is the final AEGIS
annotation. The protein FASTAs contain all translated proteins and the main
representative protein set.

## Core Evidence Outputs

```text
${output_dir}/
  assembly_masked.EDTA.fasta
  liftoff_previous_annotations.gff3
  unmapped_features.txt
  augustus.hints.gff3
  genemark.gtf
  genemark_supported.gtf
  braker.gff3
```

`assembly_masked.EDTA.fasta` is the EDTA-masked assembly. Liftoff outputs are
the transferred previous annotation and the list of unmapped features. BRAKER3
publishes AUGUSTUS and GeneMark evidence at the top level for compatibility and
for direct inspection.

Merged transcript evidence is published at the top level:

```text
${output_dir}/
  merged_hisat2_stringtie_*.gtf
  merged_star_stringtie_*.gtf
  merged_star_psiclass_*.gtf
  merged_minimap2_stringtie_long_reads_*.gtf
```

These files summarize short-read STAR/HISAT2/PsiCLASS/StringTie evidence and
long-read Minimap2/StringTie evidence.

## EGAPx Outputs

```text
${output_dir}/egapx/
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
TITAN filenames. `egapx_out/` is the nested EGAPx output directory.

## Additional Annotations

```text
${output_dir}/additional_annotations/
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

## Optional Final Annotation Source

```text
${output_dir}/final_annotations/mikado/
  final_mikado_annotation.gff3
  mikado.loci.gff3
  mikado.subloci.gff3
  versions.yml
  intermediate/
  transdecoder/
```

Mikado is an optional alternative final annotation source. The AEGIS annotation
remains under `aegis_outputs/`.

## Functional Annotation

```text
${output_dir}/Diamond2GO_outputs/
  final_annotation_proteins_all.diamond2go.tsv
  final_annotation_proteins_main.diamond2go.tsv
  versions.yml
${output_dir}/EggNOG_outputs/
  final_annotation_proteins_all.emapper.annotations
  final_annotation_proteins_main.emapper.annotations
  final_annotation_proteins_all.emapper.seed_orthologs
  final_annotation_proteins_main.emapper.seed_orthologs
  final_annotation_proteins_all.emapper.orthologs
  final_annotation_proteins_main.emapper.orthologs
  final_annotation_proteins_all.emapper.annotations.xlsx
  final_annotation_proteins_main.emapper.annotations.xlsx
  versions.yml
${output_dir}/InterProScan_outputs/
  final_annotation_proteins_all.tsv
  final_annotation_proteins_main.tsv
  final_annotation_proteins_all.gff3
  final_annotation_proteins_main.gff3
  final_annotation_proteins_all.json
  final_annotation_proteins_main.json
  versions.yml
```

Diamond2GO is part of the default graph. eggNOG-mapper and InterProScan are
optional and only appear when enabled.

## Quality, Validation And Provenance

```text
${output_dir}/quality_report/
  titan_multiqc_report.html
  titan_multiqc_report_data/
  agat_stats/
  busco/
  expression_validation/
  final_annotation_sources/
  ncrna_annotations/
  omark/
  sqanti3/
  versions.yml
${output_dir}/validation/
  final_annotation_validation.json
  final_annotation_validation.txt
  versions.yml
${output_dir}/provenance/
  evidence_manifest.json
  additional_annotations_manifest.json
  versions.yml
```

`quality_report/titan_multiqc_report.html` is the main QC entry point. The
validation directory reports structural checks on the final annotation.
Provenance manifests record inputs, selected outputs, versions and checksums.

## Intermediate And Runtime Outputs

```text
${output_dir}/intermediate_files/
  aegis/
  braker3/
  evidence_data/
    EDTA/
    RNAseq_alignments/
      HISAT2/
      STAR/
      minimap2/
    RNAseq_data/
      trimmed_data/
    transcriptomes/
      STAR_PsiCLASS/
      StringTie/
    hisat2_databases/
    star_databases/
  liftoff/
  salmon_strand/
${output_dir}/tmp/
${output_dir}/nextflow_reports/
  <run_name>.dag.html
  progress.log
```

`intermediate_files/` and `tmp/` are controlled by `--publish_intermediates`.
They are useful for debugging, manual inspection and downstream reuse, but the
main deliverables are `aegis_outputs/`, `quality_report/`, `validation/` and
`provenance/`.

The Nextflow `work/` directory and `.nextflow.log` are runtime state, not
published biological outputs. Keep them for resume/debugging, but do not treat
them as deliverables.
