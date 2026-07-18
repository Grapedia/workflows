# TITAN documentation cleanup plan

Goal: turn TITAN's documentation into something safe to hand to external
developer colleagues and users, and to cite/link from a scientific
publication — no development-session traces (dated debugging notes, task IDs,
AI-agent prompts, audit logs), no duplicated/contradictory content, a single
place to start reading, and one diagram that actually shows every step,
every tool and every option currently in the graph.

This file is a plan, not the rewrite itself. Completed items are marked
`[done]`. Each section names the exact file(s), what is wrong with them
today, and what to do. Priority order matches the numbered sections.

---

## 0. Already done

- [done] `to_do_add.md` deleted (`git rm`). It was a personal development
  roadmap/status log (13 phases, dated "Statut TITAN codex-dev" entries,
  milestone tracking) with zero value to an external reader.
- [done] Section 1 development-session artifacts removed. Durable architecture
  rationale was rewritten into `docs/development/ARCHITECTURE.md`.
- [done] Section 2 broken/stale README references fixed.
- [done] Section 3 developer references cleaned and linked from
  `CONTRIBUTING.md`.

---

## 1. Delete — pure development-session artifacts

These files exist only because an AI coding agent (Codex/Claude) needed
operating instructions or produced a session transcript. They reference a
throwaway dev VM path (`/home/vmadmin/...`), a branch name from a past
session (`codex/titan-hardening`), specific dated incidents, and internal
task IDs (`P0-00x`, `P1-00x`). None of it is meaningful to a developer
colleague reading the repo fresh, and it actively signals "this was
AI-generated" to a publication reviewer.

| File | What it is | Action |
| --- | --- | --- |
| `dev.md` | Copy-paste bootstrap prompt for a Codex worktree session (`.claude-worktrees/codex-dev`), references a specific past incident | **Delete** |
| `prompt.md` | Full AI-agent operating instructions ("Prompt de developpement pour TITAN"), references `/home/vmadmin/atcg-rnaseq` | **Delete** |
| `docs/development/audit.md` | Dated (2026-07-16) session log: git commands run, directory listing, "P1-XXX" task references, branch `codex/titan-hardening` | **Delete** |
| `docs/development/p0-hardening.md` | Dated session summary of "P0" tasks, references the two files above | **Delete** |

Before deleting `docs/development/audit.md` and `p0-hardening.md`, check
whether any *durable* fact only lives there (e.g. a design decision with no
other record). Section 3 below folds the few facts worth keeping into
`docs/development/architecture-audit.md` (see next section) so nothing
technical is lost, just the session narrative around it.

Also worth a decision: `docs/development/architecture-audit.md` itself
(282 lines) is written the same way (dated, "P1-002 did X", "Highest-risk
findings"), but roughly a third of it is durable architectural rationale
(why `main.nf` is a thin entrypoint, why EDTA/EGAPx are mandatory, why
long-read detection comes from the samplesheet and not a flag). Recommend
**rewriting it as `docs/development/ARCHITECTURE.md`**: present tense, no
dates, no P1-XXX numbering, organized by topic ("Why main.nf stays thin",
"Why EDTA and EGAPx are mandatory", "Evidence channel contracts", "Container
pinning policy") instead of as an audit trail. Everything session-flavored
(git log excerpts, "Executive summary" framed as a review of past work)
gets cut.

---

## 2. Delete or fix — broken/stale references

- [done] **`README.md` "Workflow Diagram" section** (`![Workflow Diagram](data_example/TITAN_diagram.jpg)`):
  this file does not exist on disk — it's a dead image link right now. Being
  replaced by the Mermaid diagram in section 4 below, so remove the `<img>`
  reference entirely once that lands rather than trying to regenerate the jpg.
- [done] **`README.md` line 7**: "Development audits and implementation notes are
  under `docs/development`" — update this pointer once section 1/3 lands
  (fewer, cleaner files there), and reword it as a contributor pointer, not
  a "here's our dev history" pointer.

---

## 3. Keep, but move/clean — genuinely useful developer reference

These have real, timeless content (not session narrative) and should stay,
just relocated/tidied so `docs/development/` reads as "how to contribute"
rather than "how the AI agent worked":

| File | Verdict |
| --- | --- |
| `docs/development/nextflow-dsl2-conventions.md` | [done] Kept as-is and linked from `CONTRIBUTING.md`. |
| `docs/development/container-locks.md` | [done] Kept with stable filename and linked from `CONTRIBUTING.md`. |
| `docs/development/architecture-audit.md` | [done] Rewritten into `docs/development/ARCHITECTURE.md` as durable rationale only. |

Recommended end state for `docs/development/`:

```
docs/development/
  ARCHITECTURE.md            (rewritten from architecture-audit.md)
  nextflow-dsl2-conventions.md
  container-locks.md
```

[done] `CONTRIBUTING.md` is now the contributor entrypoint and links to the
three stable `docs/development/` references instead of duplicating them.

---

## 4. `README.md` — restructure

Current state after sections 1-3: `README.md` is cleaner, but still has one
flat user-facing file mixing quick start, detailed per-tool reference
sections, a full outputs table, and troubleshooting. There is still no
pipeline diagram near the intro, so a reader has no compact way to see what
the pipeline actually does before scrolling through setup instructions.

Plan:

1. [done] **Move "Developer Quality Contract" to `CONTRIBUTING.md`**
   (see section 3). A user citing TITAN in a methods section does not need
   to know the module-authoring conventions.
2. **Add the pipeline diagram right after the intro** (section 5 below),
   before "Quick Start" — a reader (colleague or reviewer) should see the
   whole graph before reading setup instructions.
3. **Move the 14 per-tool paragraphs (lines 224-305) to a new
   `docs/reference/tools.md`**, one section per tool, keep the same content
   (it's good, accurate, current — just too long for the README itself).
   Replace them in `README.md` with a compact table: tool name, flag,
   default, one-line purpose, link to the detail section. This is the
   biggest single length reduction available (roughly 80 lines → a
   15-row table).
4. **Keep in `README.md`**: intro/contributors, Current Contract, Quick
   Start, Requirements, Input Files, RNA-seq/Protein samplesheet formats,
   EGAPx input, Profiles, the compact tool table from point 3, the Outputs
   table (it's the single most useful reference for a new user and is
   already reasonably tight), Resume/Re-runs, Troubleshooting, Limitations,
   Tool References (citations — useful for a methods section).
5. **Move "Validation and CI" section to `CONTRIBUTING.md`** — it's about
   running the dev test suite, not about using the pipeline.

Target: `README.md` under ~250 lines, `docs/reference/tools.md` holding the
detailed per-tool behavior, `CONTRIBUTING.md` holding the module-authoring
contract + test/validation instructions.

---

## 5. New: full pipeline diagram

This is the "graphe total de toutes les étapes avec outils, options, comment
c'est liés" the README needs. Use a Mermaid flowchart (renders natively on
GitHub, stays text-diffable in git — no binary image to regenerate and go
stale like the dead `TITAN_diagram.jpg`).

Draft below, organized by subgraph to stay readable. Optional
branches (anything gated by a `--run_*` flag defaulting to `false`, plus
`run_helixer`/`run_eggnog_mapper`/`run_interproscan`/`run_busco`/
`run_omark` which also default `false`, and `run_transdecoder` which
defaults `true` but only matters when `run_mikado true`) are styled
distinctly from the mandatory path. This draft should be reviewed against
`workflows/titan.nf` / `subworkflows/*.nf` once more before it goes into the
README — it was built from today's audit of the graph, not re-verified line
by line against the current DSL.

```mermaid
flowchart TD
    classDef mandatory fill:#dbe9ff,stroke:#3366cc,color:#111
    classDef optional fill:#fff3cd,stroke:#cc9900,color:#111,stroke-dasharray: 4 3
    classDef qc fill:#e6f4ea,stroke:#2e7d32,color:#111

    NEW[new_assembly FASTA]:::mandatory
    PREV[previous_assembly + previous_annotations]:::mandatory
    RNASEQ[RNAseq_samplesheet]:::mandatory
    PROT[protein_samplesheet]:::mandatory
    EGAPXIN[egapx_paramfile]:::mandatory

    subgraph EV["Evidence generation (mandatory)"]
        LIFTOFF[Liftoff\nprevious annotation to new assembly]:::mandatory
        EDTA[EDTA\nrepeat masking]:::mandatory
        EGAPX[EGAPx\nNCBI annotation pipeline]:::mandatory
        FASTP[fastp trimming]:::mandatory
        STARIDX[STAR / HISAT2 / Minimap2\ngenome indices]:::mandatory
        STARALN[STAR + StringTie]:::mandatory
        STARPSI[STAR + PsiCLASS]:::mandatory
        HISAT[HISAT2 + StringTie]:::mandatory
        MM2[Minimap2 + StringTie\nlong reads, if present]:::mandatory
        BRAKER[BRAKER3\nAUGUSTUS + GeneMark]:::mandatory
        SALMON_STRAND[Salmon\nstrand inference]:::mandatory
    end

    subgraph NCRNA["ncRNA branches (optional, parallel to evidence generation)"]
        TRNA["tRNAscan-SE\n--run_trnascan"]:::optional
        RFAM["Infernal / Rfam\n--run_rfam --rfam_data_dir"]:::optional
        HELIXER["Helixer ab initio\n--run_helixer --helixer_model_dir"]:::optional
        FLAIR["FLAIR long-read isoforms\n--run_flair"]:::optional
    end

    NEW --> LIFTOFF
    PREV --> LIFTOFF
    NEW --> EDTA
    EGAPXIN --> EGAPX
    RNASEQ --> FASTP --> STARIDX --> STARALN & STARPSI & HISAT & MM2
    PROT --> BRAKER
    EDTA --> BRAKER
    STARALN & HISAT --> SALMON_STRAND
    NEW --> TRNA
    NEW --> RFAM
    EDTA --> HELIXER
    MM2 --> FLAIR
    LIFTOFF -.splice-junction correction.-> FLAIR

    subgraph AEGIS_MERGE["AEGIS integration (mandatory)"]
        AEGIS["AEGIS merge\nrenames to Vitvi... scheme"]:::mandatory
        FINALGFF[final_annotation.gff3\nfinal_annotation_proteins_*.fasta]:::mandatory
    end

    LIFTOFF --> AEGIS
    EGAPX --> AEGIS
    BRAKER --> AEGIS
    STARALN & STARPSI & HISAT & MM2 --> AEGIS
    HELIXER -.optional evidence.-> AEGIS
    FLAIR -.optional evidence.-> AEGIS
    AEGIS --> FINALGFF

    subgraph FUNC["Functional annotation (AEGIS proteins)"]
        DIAMOND[Diamond2GO]:::mandatory
        EGGNOG["eggNOG-mapper\n--run_eggnog_mapper"]:::optional
        IPS["InterProScan\n--run_interproscan"]:::optional
    end

    FINALGFF --> DIAMOND
    FINALGFF -.-> EGGNOG
    FINALGFF -.-> IPS

    subgraph LNC["lncRNA candidates (optional)"]
        CPAT["CPAT-plant filter\n--run_lncrna --cpat_model_dir"]:::optional
    end

    FINALGFF --> CPAT
    TRNA -.exclude tRNA overlap.-> CPAT
    RFAM -.exclude ncRNA overlap.-> CPAT
    STARALN & STARPSI & HISAT & MM2 -.candidate transcripts.-> CPAT

    subgraph MIKADO_BRANCH["Mikado alternative final source (optional)"]
        MIKPREP["mikado configure + prepare\n--run_mikado"]:::optional
        TD["TransDecoder\n--run_transdecoder (default true, needs run_mikado)"]:::optional
        MIKSER["mikado serialise --orfs"]:::optional
        MIKPICK["mikado pick"]:::optional
        MIKGFF["final_mikado_annotation.gff3"]:::optional
    end

    LIFTOFF & EGAPX & BRAKER & HELIXER & FLAIR --> MIKPREP
    STARALN & STARPSI & HISAT & MM2 --> MIKPREP
    MIKPREP --> TD --> MIKSER --> MIKPICK --> MIKGFF

    subgraph QUALITY["Quality report (mandatory shell, contents mostly optional)"]
        BUSCO["BUSCO\n--run_busco (Vitis default: on)"]:::optional
        OMARK["OMArk\n--run_omark --omark_data_dir"]:::optional
        AGATSTATS[AGAT structural stats]:::mandatory
        NCRNAQC[ncRNA count summary]:::mandatory
        SQANTI["SQANTI3\n--run_sqanti3"]:::optional
        EXPR["Expression support validation\nSalmon on final transcripts (default: on)"]:::mandatory
        SRCQC["AEGIS vs Mikado comparison"]:::optional
        MULTIQC[MultiQC HTML report]:::mandatory
    end

    FINALGFF --> BUSCO & AGATSTATS & EXPR
    DIAMOND & EGGNOG & IPS --> MULTIQC
    TRNA & RFAM --> NCRNAQC
    MM2 & FLAIR --> SQANTI
    FINALGFF --> SQANTI
    FINALGFF --> OMARK
    FINALGFF & MIKGFF --> SRCQC
    BUSCO & OMARK & AGATSTATS & NCRNAQC & SQANTI & EXPR & SRCQC --> MULTIQC
```

Legend to include next to the diagram in the README:

* Blue, solid border = always runs.
* Yellow, dashed border = optional, off by default unless noted, enabled
  with the `--run_*` flag shown.
* Green = quality-report step.

This diagram intentionally omits fine-grained sub-steps (e.g. individual
STAR/HISAT2 per-sample fan-out, the four StringTie stranded/unstranded
variants) to stay readable — it documents the tool-level graph, not the
per-task DAG (Nextflow's own `-with-dag` output already covers that at
full resolution per run).

---

## 6. `docs/user/installation.md` — gap to fill

589 lines, well-organized (14 numbered sections), no development-session
traces found. Genuinely good user documentation. **One real gap**: sections
8-10 cover eggNOG-mapper/Helixer/InterProScan offline-data setup, but there
is no equivalent section for the 8 tools added since (tRNAscan-SE, Rfam,
lncRNA/CPAT, Mikado, TransDecoder, FLAIR, SQANTI3, OMArk) — a new user
following this guide today would have no instructions for
`scripts/download_rfam_data.sh` / `scripts/download_omark_data.sh` or the
corresponding `--prepare-rfam-data` / `--prepare-omark-data` launcher flags,
even though `README.md` documents the tools themselves.

Action: add sections 8a-8h (or renumber) following the exact pattern of the
existing eggNOG-mapper/Helixer/InterProScan sections (what the tool does,
how to fetch offline data once, which flags enable it, link to the
`README.md` tool-reference entry once section 4 lands).

---

## 7. `data/slurm_apptainer.config` — trim dev-session narrative from comments

This file is tracked in git (confirmed during the codex-dev merge work —
it is **not** covered by `.gitignore` despite living under `TITAN/data/`,
because it was added before that ignore pattern existed). It will ship to
every colleague who clones the repo, so its comments are effectively public
documentation of the Colmar cluster deployment, not private scratch notes.

Most of it is legitimate, valuable operational tuning rationale (why
`process_trim` gets less memory than `process_alignment`, why some labels
need a narrower `--nodelist`). Keep the *facts*. Rewrite the *framing*:
comments like "used to share process_alignment and that alone was enough to
starve the Slurm queue for hours on job 532274 (2026-07-17)" should become
"fastp only needs a couple of threads and a few GB of RAM, unlike real
alignment/index jobs — giving it `process_alignment`'s 96GB budget wastes
node capacity it doesn't need." Same fact, no incident-report framing, no
job ID, no date.

Same treatment for the one comment in `conf/slurm.config` referencing
"confirmed empirically" without further detail — keep the conclusion, drop
the "confirmed empirically" framing (either state it as fact, or footnote
briefly why, without narrating the debugging process).

This is a low-risk, mechanical pass: read every comment in both files,
keep the technical conclusion, cut anything written like a debugging
session note (dates, job IDs, "we found that...", "confirmed empirically").

---

## Suggested order of execution

1. Section 1 deletions (`dev.md`, `prompt.md`, `docs/development/audit.md`,
   `docs/development/p0-hardening.md`) — fast, unambiguous, zero risk once
   confirmed.
2. Section 7 (`data/slurm_apptainer.config` / `conf/slurm.config` comment
   cleanup) — mechanical, no structural change, safe to do independently.
3. Section 3 (`architecture-audit.md` → `ARCHITECTURE.md` rewrite,
   `CONTRIBUTING.md` creation).
4. Section 4 + 5 together (`README.md` restructure needs the diagram to
   land in the same pass, otherwise the README is briefly missing the thing
   it's being restructured to showcase).
5. Section 6 (`docs/user/installation.md` new tool sections) — independent,
   can happen any time after section 4 gives it something to link to.

Everything above is additive/reversible through git except the section 1
deletions, which is why those are listed as a plan here rather than already
executed alongside `to_do_add.md`.
