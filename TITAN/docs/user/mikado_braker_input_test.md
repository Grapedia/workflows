# Test: feeding Mikado the TSEBRA-filtered BRAKER3 tracks

Follow-up to [unsupported_genes.md](unsupported_genes.md). That analysis found that
82 % of the 2,478 genes proposed to be set aside are raw BRAKER3 models, and that
`mikado_prepare` receives the **raw** `augustus.hints.gff3` and `genemark.gtf` instead
of the TSEBRA-selected `braker.gff3` and GeneMark's `genemark_supported.gtf`. This
test measures what changes if Mikado receives the filtered tracks.

Status: **run 2026-09-29, verdict: do not switch production as is** (see [Results](#results)).
Production (`titan.nf`, modules) is unchanged.

## Design

`mikado_input_test.nf` (project root, so `projectDir` and `nextflow.config` are the
production ones) re-runs only

`mikado_prepare` -> `transdecoder_longorfs` -> `transdecoder_predict` -> `mikado_serialise` -> `mikado_pick`

with the production modules and containers, on the evidence tracks that a finished
run has already published (`data/titan_prod_out/04_evidence`, `03_additional_annotations`).
It then extracts the longest-isoform proteins from the Mikado loci and runs the
production `busco` and `agat_stats` modules on them.

| Arm | Track 1 (`braker_augustus` slot) | Track 2 (`braker_genemark` slot) |
|---|---|---|
| `raw` (control) | `augustus.hints.gff3` | `genemark.gtf` |
| `tsebra` (test) | `braker.gff3` | `genemark_supported.gtf` |

Everything else is identical, including labels and priority scores, so the only
variable is the content of the two BRAKER3 tracks. Labels are deliberately not renamed
in the test (`braker.gff3` already contains GeneMark-derived models); rename them if
this becomes the production behaviour.

The **control arm is essential**. It must reproduce the production Mikado output;
otherwise the harness differs from production and the `tsebra` arm cannot be
interpreted. (The production work directory was cleaned, so the published files are
the only source.) Known small differences from production: the four HISAT2 slots
get empty placeholders as in production (`run_hisat2 = false`), and the flair track is
the published, empty `flair_isoforms.gtf`.

Not tested here: AEGIS renaming/tidying, ncRNA/lncRNA calling and the downstream
functional annotation. The comparison is at the level of the Mikado loci. Also,
`genemark_supported.gtf` keeps every `gmst` model (they are RNA-seq-derived, so all are
"supported") and only filters the `GeneMark.hmm3` ones (492,372 -> 272,756 lines).

## Run

```bash
sbatch scripts/mikado_input_test/launch_arm.sh raw
sbatch scripts/mikado_input_test/launch_arm.sh tsebra
```

Each arm submits its own Slurm tasks (Mikado prepare/serialise/pick, TransDecoder,
BUSCO, AGAT) with the resources of `data/slurm_apptainer.config`
(`process_transcriptome`: 16 CPUs, 96 GB; BUSCO: `process_aegis` 16 CPUs, 96 GB).
Outputs: `data/mikado_input_test/<arm>/` (same layout as a TITAN output dir; Nextflow
work in `data/mikado_input_test/work_<arm>`). Re-running the same command resumes.

Nextflow needs Java 11+; the launcher loads `Java/17.0.13`, `nextflow/24.04.3`,
`apptainer` and `python` modules like the production launcher. (In an
non-interactive shell `module load` may not take effect; use
`eval "$($LMOD_CMD bash load ...)"`.)

Compare when both arms are finished:

```bash
python3 scripts/mikado_input_test/compare_arms.py \
  --raw-dir data/mikado_input_test/raw --tsebra-dir data/mikado_input_test/tsebra \
  --prod-mikado data/titan_prod_out/04_evidence/mikado/final_mikado_annotation.gff3 \
  --prod-gff3 data/titan_prod_out/01_final_annotation/primary/final_annotation.gff3 \
  --classes data/unsupported_genes_analysis/classes/gene_classes.tsv \
  --out-dir data/mikado_input_test/comparison
```

The script was checked on a control (production Mikado output given as both arms):
every production gene is retained and no new gene appears.

## Decision criteria (fixed before running)

These are a proposal; change them before running, not after.

1. **The control is valid**: the `raw` arm reproduces >= 99 % of the production Mikado
   loci (`loci_identical_to_production_mikado`) and its BUSCO is within 0.3 points of
   the production Mikado proteins.
2. **The test removes the right genes**, using the classes of the unsupported-genes
   analysis on the production genes:
   - loses >= 50 % of the 2,478 `no_signal` + `te_overlap_unsupported` genes;
   - loses <= 1 % of the `expressed` genes and <= 5 BUSCO genes;
   - loses <= 10 % of the `conserved_or_functional` genes (they are probably real).
3. **Completeness does not degrade**: BUSCO complete drops by <= 0.3 points relative
   to the control arm.
4. **Few new genes**: genes present in the `tsebra` arm but overlapping no production
   gene are <= 2 % of the total (new calls have no evidence table yet).

"Clearly better" = 1 to 4 all met. If 2 is met but 3 fails, the raw models are buying
completeness and the answer is a compromise (for example keep the raw tracks with a
lower priority score), not a switch. If the `raw` arm fails 1, stop and fix the harness.


## Results

Both arms ran with 24-CPU Mikado tasks (`resources.config`). Wall time per arm: about 4 h,
of which TransDecoder Predict alone is 2 h 17 min (`raw`) and 2 h 22 min (`tsebra`); it is
single-threaded (one Perl `score_CDS_likelihood_all_6_frames.pl` process), so the 24 CPUs
only help `mikado_serialise` and `mikado_pick`. Full tables: `data/mikado_input_test/comparison/`.

| | production Mikado | `raw` (control) | `tsebra` (test) |
|---|---:|---:|---:|
| Genes | 55,766 | 55,767 | **49,700** |
| mRNAs | 57,542 | 57,543 | 51,528 |
| Loci identical to production | - | 55,323 (99.2 %) | 44,248 |
| Genes overlapping no production gene | - | 400 | 1,151 |
| Single-exon genes | - | 15,641 | 14,584 |
| BUSCO complete (eudicotyledons_odb12.2, n=1,990) | 99.8 % (final, after AEGIS) | 99.8 % | 99.9 % |

Fate of the production protein-coding genes (a gene is "retained" if >= 50 % of its CDS
overlaps a CDS of the arm, same strand), by class of
[unsupported_genes.md](unsupported_genes.md):

| Class | Genes | Lost in `raw` (control) | Lost in `tsebra` |
|---|---:|---:|---:|
| expressed | 41,290 | 392 (0.9 %) | **4,023 (9.7 %)** |
| conserved_or_functional | 10,308 | 0 | **2,429 (23.6 %)** |
| te_protein | 1,140 | 0 | 411 (36.1 %) |
| te_overlap_unsupported | 1,215 | 11 (0.9 %) | 709 (58.4 %) |
| liftoff_only | 539 | 0 | 30 (5.6 %) |
| junction_only | 11 | 0 | 2 |
| no_signal | 1,263 | 0 | 781 (61.8 %) |
| BUSCO genes (of 2,022) | | 0 | **0** |

### Verdict against the criteria fixed beforehand

| Criterion | Result | |
|---|---|---|
| 1. Control reproduces production (>= 99 % identical loci, BUSCO within 0.3) | 99.2 %, BUSCO 99.8 % | met |
| 2a. Loses >= 50 % of the 2,478 set-aside genes | 1,490 lost (60.1 %) | met |
| 2b. Loses <= 1 % of expressed genes and <= 5 BUSCO genes | 9.7 % expressed, 0 BUSCO | **not met** |
| 2c. Loses <= 10 % of conserved_or_functional | 23.6 % | **not met** |
| 3. BUSCO does not drop by more than 0.3 | 99.9 % vs 99.8 % | met |
| 4. New genes <= 2 % of the total | 1,151 (2.3 %) | **not met** (barely; the control alone gives 0.7 %) |

**The test is not "clearly better", so production should not be switched as is.** It is a
filter, and a blunt one: it removes 60 % of the targeted models but also a large part of
the BRAKER models that carry expression or function.

### What the losses look like

- **They are almost all BRAKER models**: 7,822 of 18,429 BRAKER-origin genes are lost
  (42 %); Helixer 0.4 %, EGAPx 0 %, Liftoff 0.2 %, RNA-seq assemblies 8.9 %.
- **AUGUSTUS drives it**: 55 % of the models whose winning transcript is a raw AUGUSTUS
  model are lost, against 27 % of the GeneMark ones. The `braker.gff3` swap matters more
  than the `genemark_supported.gtf` swap; the two changes were not separated.
- **The selection points the right way**: within BRAKER-origin genes of the same class,
  the lost ones are weaker than the kept ones. Expressed: median protein 147 vs 298 aa,
  intron support 3 % vs 26 %, non-*Vitis* homolog 32 % vs 59 %, in PN40024 5.1 19 % vs 45 %,
  median max TPM 1.9 vs 3.1. Conserved: protein 142 vs 297 aa, OMAMER 44 % vs 72 %.
- **But not all of them are junk**: among the 4,023 lost expressed genes, 26 % have a max
  TPM >= 5 and 64 % have a functional hit; 306 combine TPM >= 5, function and a non-*Vitis*
  homolog. Among the 2,429 lost conserved genes, 583 have function, OMAMER family and a
  homolog. About 890 genes look real (upper bound of a crude definition).
- **They are truly gone, not replaced**: for the lost expressed genes, only 30 % have any
  same-strand arm gene at the locus (50 % on either strand).
- **BUSCO cannot see it**: 99.9 % with about 6,000 fewer genes. BUSCO is insensitive to
  this kind of loss and does not validate either option.

Reading: many of the "expressed" BRAKER-only genes probably owe their TPM to reads that
are not from a distinct gene (overlapping neighbours, unspliced or antisense reads:
only 3 % of the lost expressed multi-exon genes have any RNA-seq intron). This is a
hypothesis; it was not tested.

### Consequences

- The classification of [unsupported_genes.md](unsupported_genes.md) is **more precise
  than the pipeline change**: it flags 2,478 genes with 0 BUSCO genes and 6-8 % overlap
  with PN40024 5.1, and it does not touch the expressed or conserved genes.
- Recommended: keep the production Mikado inputs and apply the classification as a flag
  (or a filtered candidate next to `high_confidence_monoexonic`).
- If a pipeline change is still wanted, test a third arm that separates the effects
  (`braker.gff3` only, or `genemark_supported.gtf` only), or keep the raw tracks with a
  lower priority score; do not adopt `tsebra` wholesale.
- Caveats: the control differs from production by 0.8 % of loci (443), so every loss
  figure includes about 1 % of baseline noise; the comparison is at the Mikado level
  (before AEGIS, ncRNA calling and functional annotation); "lost" uses a 50 % CDS
  overlap rule, so a locus that changes model but keeps the CDS is retained.

## If a change is adopted: changing production

The switch is in `workflows/titan.nf` (the two `evidence_data.braker_*` arguments of
`mikado_prepare`) and `modules/mikado.nf` (labels `braker_augustus`/`braker_genemark`
and their scores 18/17). Steps: make the choice a parameter (default = current
behaviour), rename the source labels, update `docs/reference/tools.md` and
`docs/user/inputs_outputs.md`, re-run the pipeline from `mikado_prepare` onward
(AEGIS, expression validation, functional annotation, BUSCO/OMArk and the audit all
depend on the gene set), then regenerate the unsupported-genes analysis on the new
run. Do not edit the production output directory in place.
