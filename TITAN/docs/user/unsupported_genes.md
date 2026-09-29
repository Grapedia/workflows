# Genes without support: how to study them and when to set them aside

Analysis of the PN40024 T2T production run (`data/titan_prod_out`, 56,744 genes),
written 2026-09-29. It answers three questions:

1. What are the ~15,000 genes with no RNA-seq expression?
2. Which of them can be called annotation errors, and why?
3. What in the pipeline produces them?

Everything below is reproducible with `scripts/run_unsupported_genes_analysis.sh`
(see [Reproduce](#reproduce)). Outputs live in `data/unsupported_genes_analysis/`.

A French HTML report that puts this analysis in the context of the overall annotation quality is
in [`docs/reports/annotation_quality_and_gene_count.html`](../reports/annotation_quality_and_gene_count.html)
(built by `scripts/build_annotation_report.py`). The optional tool that marks and filters a GFF3 with
these classes is `scripts/flag_unsupported_genes.py` (see [Marking and filtering](#marking-and-filtering-a-gff3)).

## Contents

- [Summary](#summary)
- [Definitions](#definitions)
- [Evidence collected per gene](#evidence-collected-per-gene)
- [Result: the 15,454 genes split into distinct populations](#result)
- [Why the flagged genes are errors](#why-the-flagged-genes-are-errors)
- [Where they come from in the pipeline](#where-they-come-from-in-the-pipeline)
- [Hypotheses that were tested and rejected](#hypotheses-that-were-tested-and-rejected)
- [The unexpressed genes that are probably real](#the-unexpressed-genes-that-are-probably-real)
- [Effect on the gene count](#effect-on-the-gene-count)
- [Recommended policy](#recommended-policy)
- [Limitations](#limitations)
- [Marking and filtering a GFF3](#marking-and-filtering-a-gff3)
- [Reproduce](#reproduce)

## Summary

- "No expression" is **15,454 genes**, but that number is not a homogeneous pool.
  978 of them are ncRNA genes that were never quantified (their TPM is exactly 0
  in every library because Salmon ran on the protein-coding transcripts), so the
  real figure is **14,476 protein-coding genes** out of 55,766.
- Only **2,478 of those (4.4 % of the protein-coding genes)** have nothing else
  behind them and can be defended as errors:
  - **1,263 `no_signal`**: no expression, no splice junction in any RNA-seq
    assembly, no functional hit, no OMAMER family, no homolog in non-*Vitis*
    eudicots, not a v4.3 gene, not on a transposable element.
  - **1,215 `te_overlap_unsupported`**: at least 50 % of the CDS lies inside an
    EDTA transposable-element annotation and, again, nothing else supports it.
- **10,308 unexpressed genes are probably real** (functional hit, OMAMER
  placement or a homolog in another eudicot). They are simply not expressed in the
  47 libraries, or their reads cannot be assigned (identical copies, mostly on
  chr00). They should be kept.
- **1,140 unexpressed genes encode real transposable-element proteins**
  (reverse transcriptase, gag-pol, integrase, transposase). They are not errors of
  prediction but they are not host genes either. Another 1,289 expressed genes
  are in the same situation. They belong in a TE track.
- The 2,478 flagged genes are **82 % raw BRAKER (AUGUSTUS/GeneMark) models**.
  Mikado is fed the *raw* AUGUSTUS and GeneMark predictions, not the TSEBRA-filtered
  set, and keeps any locus that nothing else competes with.
- Setting the 2,478 aside leaves about **53,300 genes**. That does **not** close the
  gap with the ~35,000 genes expected for grapevine; it is not the main explanation.

## Definitions

| Term | Meaning |
|---|---|
| expressed | TPM >= 0.5 in at least one of the 47 Salmon libraries (same rule as the TITAN audit) |
| function | any hit in Diamond2GO, eggNOG-mapper or InterProScan, or a MapMan (Mercator4) bin other than 99 "not assigned" |
| OMAMER hit | protein placed in a hierarchical orthologous group by OMArk's OMAMER search |
| non-*Vitis* homolog | DIAMOND hit (e < 1e-5, >= 35 % identity, >= 50 % query coverage) against UniProt eudicot proteins **excluding all *Vitis* entries** |
| junction support | at least one intron of the gene is reproduced exactly by an intron of a StringTie/PsiCLASS (STAR) or long-read assembly |
| TE overlap | >= 50 % of the CDS overlaps an EDTA TE annotation |
| TE protein hit | a functional hit whose text names a TE protein (transposase, reverse transcriptase, gag-pol, integrase, copia/gypsy, MULE, hAT ...) |

A non-*Vitis* homolog is used on purpose: a hit against another *Vitis* protein may
just be another *Vitis* annotation of the same (possibly wrong) model.

## Evidence collected per gene

`data/unsupported_genes_analysis/gene_evidence_table.tsv`, one row per gene, built
by three scripts:

| Script | Adds |
|---|---|
| `scripts/build_gene_evidence_table.py` | gene structure (exons, CDS, UTR), Mikado provenance (`alias` of the winning model, joined by coordinates: 100 % matched), TPM summary, functional hits (Diamond2GO, eggNOG, InterProScan, MapMan), OMAMER family, Liftoff carry-over, protein composition, EDTA overlap, neighbourhood (overlaps, distances), TE-protein hits |
| `scripts/add_intron_support.py` | how many introns are found in RNA-seq assemblies (independent of Salmon) |
| `scripts/add_homology_features.py` | TSEBRA overlap, DIAMOND hits (self, Vitales+Swiss-Prot, non-*Vitis* eudicots), paralog/fragment features |

`scripts/salmon_gene_numreads.py` adds the total Salmon `NumReads` per gene.

Sanity check of the table: the number of genes with no expression, no function and
no Liftoff ID is 3,416, which matches `ANNOTATION_QUALITY_AUDIT.md` exactly.

## Result

Protein-coding genes only (55,766). `scripts/classify_unsupported_genes.py` assigns
each gene exactly one class, in this priority order:

| Class | Genes | % | Median protein (aa) | BUSCO genes | Overlap with PN40024 5.1 |
|---|---:|---:|---:|---:|---:|
| expressed (reference) | 41,290 | 74.0 | 288 | 2,021 | 66.8 % |
| conserved_or_functional (unexpressed) | 10,308 | 18.5 | 190 | 1 | 19.5 % |
| te_protein (unexpressed) | 1,140 | 2.0 | 250 | 0 | 12.7 % |
| **te_overlap_unsupported** | **1,215** | 2.2 | 114 | 0 | 8.2 % |
| liftoff_only (v4.3 gene, nothing else) | 539 | 1.0 | 107 | 0 | 22.3 % |
| junction_only | 11 | 0.0 | 93 | 0 | 18.2 % |
| **no_signal** | **1,263** | 2.3 | 115 | 0 | 6.1 % |

Origin of the winning Mikado model (`alias`), by class:

| Class | BRAKER (raw) | Helixer | EGAPx | Liftoff v4.3 | RNA-seq assembly |
|---|---:|---:|---:|---:|---:|
| expressed | 8,434 | 4,158 | 16,813 | 6,139 | 5,746 |
| conserved_or_functional | 6,849 | 1,179 | 215 | 2,016 | 49 |
| te_protein | 1,025 | 28 | 16 | 59 | 12 |
| te_overlap_unsupported | 1,046 | 21 | 0 | 97 | 51 |
| liftoff_only | 89 | 14 | 0 | 436 | 0 |
| no_signal | 977 | 277 | 0 | 0 | 9 |

EGAPx and the RNA-seq assemblies almost never appear in the unsupported classes:
they are built from the same RNA-seq, so "expressed" is close to true by
construction for them. The expression criterion therefore only discriminates among
ab initio (BRAKER, Helixer) and Liftoff models, and that is where the problem is.

## Why the flagged genes are errors

There is no proof that a gene does not exist. The argument is that every
independent line of evidence agrees, and that the flagged set differs from genes
known to be real in every dimension, while control sets do not.

1. **Not just below a threshold.** Total Salmon `NumReads` over all 47 libraries
   (`classes/salmon_reads_by_class.tsv`): `no_signal` median **1 read**, 41 % with
   zero reads, 89 % with fewer than 10 reads in total. Expressed genes: median
   2,353 reads, 84 % with >= 100. The absence of expression is an absence of reads,
   not a TPM cutoff artefact.
2. **No splice junction.** 96 % of the unexpressed multi-exon genes (98.5 % of
   those with no function either) have no intron reproduced by any RNA-seq assembly,
   against 30 % for expressed genes. This uses different software from Salmon.
3. **No conservation.** Among the 3,192 genes with no expression, no function and no
   MapMan bin, only 1 % have a non-*Vitis* eudicot homolog and 4 % an OMAMER family,
   against 73 % and 78 % for expressed genes. A real protein-coding gene in a flowering plant nearly always
   has a relative somewhere. Those that do not are either extremely lineage-specific
   or not genes.
4. **Not a real gene of the reference set.** No BUSCO ortholog is flagged
   (0 of 2,022 BUSCO genes). Overlap with the independent PN40024 5.1 annotation is
   6.1 % for `no_signal` and 8.2 % for `te_overlap_unsupported`, against 66.8 % for
   expressed genes. 5.1 is itself RNA-seq-driven, so this is not fully independent
   of expression; it is a consistency check, not a proof.
5. **Small, ab initio-shaped models.** Median protein 115 aa (expressed: 288 aa),
   33 % shorter than 100 aa (expressed: 13 %), 72 % raw ab initio models that the
   TSEBRA step did not keep (expressed: 15 %).

For `te_overlap_unsupported` the explanation is more direct: the model sits on a
transposable element and has no host-gene evidence.

**Confidence.** High for `no_signal` and `te_overlap_unsupported` taken together.
Moderate for `liftoff_only` (539 genes: the v4.3 annotation is a prior, not evidence;
they have the same shape as `no_signal`, 93.7 % with fewer than 10 reads, but were
kept by a previous curated annotation). `junction_only` (11 genes) is too small to
decide.

## Where they come from in the pipeline

`modules/mikado.nf` feeds Mikado with `braker_augustus_gff3` and `braker_genemark_gtf`,
wired in `workflows/titan.nf` from `braker3_results.augustus_gff` and
`braker3_results.genemark_gtf`, i.e. `augustus.hints.gff3` (83,255 transcripts) and
`genemark.gtf`. These are the **raw** predictions. The TSEBRA-selected set is
`braker.gff3` (52,414 transcripts) and GeneMark's own filtered set is
`genemark_supported.gtf`; neither reaches Mikado (in `workflows/titan.nf` only the
two raw channels are passed on; `braker_gff3` and `braker_genemark_supported_gtf`
are exposed by `subworkflows/generate_evidence_data.nf` but are only published).

Measured consequence:

- 18,429 final genes have a raw BRAKER model as the winning Mikado transcript.
  Only 25 % of them overlap (>= 50 % of the CDS, same strand) a `braker.gff3`
  model. The rest were rejected or never proposed by TSEBRA.
- The winners of a Mikado locus are kept even when nothing else competes with them:
  13,697 superloci give 56,744 genes. A lone ab initio model with no RNA-seq, protein
  or junction evidence is therefore reported unchanged.
- Among unexpressed genes 23.6 % are in the TSEBRA set, among expressed 61 %.

This is an association, not a proven cause: I did not re-run Mikado with the
TSEBRA set. But it is the simplest explanation of why the unsupported models are
overwhelmingly BRAKER (82 % of the flagged set, against 20 % of all genes).

Caution before changing it: the raw models probably also contribute to the very
high BUSCO completeness (99.8 % vs 97.9 % for PN40024 5.1). Any change must be
tested by re-running Mikado and re-measuring BUSCO/OMArk.

## Hypotheses that were tested and rejected

Prevalence of each signature in `no_signal`, against expressed genes
(`classes/signature_summary.tsv`):

| Hypothesis | `no_signal` | expressed | Verdict |
|---|---:|---:|---|
| Fragment of a paralog (partial hit to a longer protein) | 10.5 % | 13.5 % | not enriched |
| Split gene: supported same-strand neighbour within 2 kb | 28.7 % | 43.8 % | *less* frequent |
| Antisense to or nested in another gene | 18.2 % | 22.7 % | not enriched |
| Same-strand CDS overlap with another gene | 6.8 % | 2.0 % | mildly enriched |
| Low-complexity protein | 3.1 % | 4.6 % | not enriched |
| chr00 (unplaced) | 6.4 % | 3.2 % | mildly enriched |
| Short ORF (< 100 aa) | 33 % | 13 % | enriched |
| Single-exon | 31 % | 28 % | not enriched |
| Raw ab initio model not kept by TSEBRA | 72 % | 15 % | strongly enriched |
| Helixer-only model | 22 % | 10 % | enriched (2x) |

So these are not mostly fragments, split genes or overlaps; they are isolated,
short, ab initio calls with no external support.

## The unexpressed genes that are probably real

10,308 genes are unexpressed but have a functional hit, an OMAMER family or a
non-*Vitis* homolog. Among them:

- 98.4 % have a functional hit, 66 % an OMAMER family, 48 % a non-*Vitis* homolog,
  35 % come from Liftoff v4.3 carry-over.
- **chr00 stands out**: 53 % of the protein-coding genes on chr00 are in this class
  (versus 17 % elsewhere). Of the 1,555 such genes on chr00, **1,376 are >= 95 %
  identical (>= 80 % coverage) to another gene**, at a median of 2.5 Mb away on chr00
  (not tandem), and 96 % originate from BRAKER. Only 25 % of their partners are
  expressed. The pattern
  fits duplicated or collapsed unplaced sequence in which reads cannot be assigned
  to a single copy; it does not fit random ORFs, because they have real functions.
  **This is a hypothesis; it was not checked at the read level** (see Limitations).
- 2,237 of the class have only weak support (function alone: no non-*Vitis* homolog,
  no OMAMER family, no v4.3 ID) and would be the first to review if the policy is
  made broader.

## Effect on the gene count

| Step | Protein-coding genes |
|---|---:|
| Final annotation | 55,766 |
| Set aside `no_signal` + `te_overlap_unsupported` (recommended) | 53,288 |
| Also remove `liftoff_only` + `junction_only` | 52,738 |
| Also remove every TE-flagged gene, expressed or not | about 41,500 |

The last line is an upper bound for the effect of TEs, not a recommendation:
7,701 expressed genes carry a TE flag, but 3,235 flagged only by a TE-like domain
(without EDTA overlap) include 13 BUSCO genes, so the domain alone is not reliable.
Among the TE flag types (unexpressed genes): 1,140 have both overlap and a TE
protein hit (91 % have an OMAMER family: genuine TE proteins), 3,174 have EDTA
overlap only (18 % OMAMER: mostly junk on TE sequence or fragments), and 1,023
have a TE domain only (80 % OMAMER: may be host genes).

The residual gap to ~35,000 is therefore not explained by unsupported calls.
Most of the extra genes are supported (expressed, or conserved). Possible remaining
causes that this analysis did not address: gene fragmentation, the ~11 % of genes in
some exon-level overlap with a neighbour, duplicated regions such as chr00, and the
35,000 figure being a conservative older baseline (see the audit, section 5.3).

## Recommended policy

1. **Do not delete, flag.** Add the class as a GFF3 attribute or ship
   `set_aside_candidates.txt` next to the primary annotation. The primary and the
   `high_confidence_monoexonic` candidates stay unchanged.
2. **Use the strict list** (`--policy strict`, 2,478 genes) as the set to set aside;
   keep `liftoff_only` for manual review.
3. **Report the 2,429 TE-protein genes separately** (1,140 unexpressed + 1,289
   expressed) as a TE track, not as errors and not as host genes.
4. **Test the pipeline cause**: re-run Mikado with `braker.gff3` and
   `genemark_supported.gtf` instead of the raw sets, and compare BUSCO, OMArk, gene
   count and the number of `no_signal` genes.
5. **Review chr00** duplicates before trusting the "unexpressed" label there.

## Limitations

- Expression is measured on 47 libraries. A gene expressed only in an untested
  tissue is unexpressed here; this is why the classes with function or
  conservation are kept.
- The `no_signal` class uses absence of evidence. It cannot detect a very
  lineage-specific real gene: expect some false positives (5.1 contains 77 of the
  1,263, 6.1 %).
- DIAMOND references (UniProt eudicots, Swiss-Prot, Vitales) can themselves contain
  ab initio proteins; a hit is evidence that a family exists, not that this model is
  correct.
- TE calls depend on EDTA and on keyword matching in functional annotations.
- Overlap with PN40024 5.1 is partly dependent on expression (5.1 is RNA-seq-based).
- "Probably an ab initio artefact" is inferred from the origin of the winning
  model and its lack of support. I did not test whether these ORFs look like random
  ORFs (codon usage, GC3), and did not re-run Mikado without the raw models.
- A direct read count from the STAR BAMs (with `bedtools multicov`) was attempted
  and abandoned: about 5 features/s on 237,352 CDS features, roughly ten hours per
  BAM. Salmon `NumReads` was used instead, which is not independent of Salmon.
  The chr00 multi-mapping hypothesis needs this read-level check (or MAPQ-filtered
  coverage) to be confirmed.
- The 978 ncRNA genes have TPM 0 because they were not quantified. TITAN's
  "27.2 % of genes with no expression" therefore includes them.

## Marking and filtering a GFF3

`scripts/flag_unsupported_genes.py` (standard library only, **not part of the pipeline**) reads
`final_annotation.gff3` and `classes/gene_classes.tsv`, never modifies its inputs and never renames a
gene. Sub-commands: `summary` (how many genes each policy removes), `mark` (adds
`titan_evidence_class=...;titan_filter_level=keep|strict|broad|te` to every `gene` line) and `filter`
(removes the genes of a policy with all their mRNA/exon/CDS/UTR children, and optionally the matching
protein FASTA records).

| Policy | Removes | Genes left (GFF3, with the 978 ncRNA) |
|---|---|---:|
| `strict` | `no_signal`, `te_overlap_unsupported` | 54,266 |
| `broad` | + `liftoff_only`, `junction_only` | 53,716 |
| `broad+te` | + `te_protein` (TE-encoded proteins, better kept in a separate TE track) | 52,576 |

`conserved_or_functional` and `expressed` genes are never removed. Checked on the production run: the
`strict` output passes `scripts/validate_final_annotation.py` (PASS, 53,288 main proteins found) with the
same warnings as the input. Unit test: `scripts/test_flag_unsupported_genes.py`. Recommendation: distribute
the **marked** annotation and let each analysis choose its policy rather than shipping a filtered file.

## Reproduce

```bash
# 1. Homology searches (DIAMOND, ~10 min on 40 threads, needs the diamond2go image)
data/unsupported_genes_analysis/diamond/run_diamond.sh all

# 2. Table, junction support, homology features, PN40024 5.1 overlap, Salmon reads, classes
scripts/run_unsupported_genes_analysis.sh \
    data/titan_prod_out data/unsupported_genes_analysis data/PN40024_5.1_on_T2T_ref.gff3
```

Python 3.6+ with pandas and bedtools are required. Files of interest:

| File | Content |
|---|---|
| `classes/gene_classes.tsv` | class, signatures and key evidence per protein-coding gene |
| `classes/class_summary.tsv` | counts per class and per origin |
| `classes/signature_summary.tsv` | signature prevalence per class |
| `classes/salmon_reads_by_class.tsv` | total Salmon reads per class |
| `classes/set_aside_candidates.txt` | the 2,478 gene IDs (strict policy) |
| `gene_evidence_table.tsv` | all 91 evidence columns for every gene |

Thresholds (`MIN_TPM 0.5`, homolog 35 % identity / 50 % coverage, TE 50 % CDS,
short ORF 100 aa) are constants at the top of the scripts.
