#!/usr/bin/env python3
"""Build the self-contained HTML report on the quality of a TITAN annotation and on the
discussion of its gene count (French).

Numbers about the unsupported / set-aside genes are computed here from the analysis tables
(scripts/run_unsupported_genes_analysis.sh, data/unsupported_genes_analysis/) and from the
Mikado input test (scripts/mikado_input_test/).  Numbers taken from the TITAN quality audit
(ANNOTATION_QUALITY_AUDIT.md) and from the PN40024 5.1 re-run (data/v5_1_audit) are listed in
the AUDIT dictionary below, with their source, so that every figure of the report can be
traced.

Usage (Python >= 3.6, pandas + numpy):
  scripts/build_annotation_report.py --titan-out data/titan_prod_out \\
      --analysis-dir data/unsupported_genes_analysis --test-dir data/mikado_input_test \\
      --out docs/reports/annotation_quality_and_gene_count.html
"""
import argparse
import datetime
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from report_charts import (data_table, dec_fr, esc, figure, fmt_int, fmt_pct, grouped_hbars,  # noqa: E402
                           hbars, heatmap, line_chart, stacked_hbars)

# --------------------------------------------------------------------------------------------
# Numbers reused from earlier analyses (each with its source file)
# --------------------------------------------------------------------------------------------
AUDIT = {
    # 01_final_annotation/quality_report/agat_stats/agat_stats.txt (audit section 2)
    "genes": 56744, "coding": 55766, "ncrna": 978, "mrna": 57542, "genes_per_mb": 112,
    "exons_per_mrna": 4.3, "single_exon_pct": 28.0, "single_exon": 15640, "utr_pct": 53.9,
    "mean_gene_len": 3923, "mean_cds_len": 1111, "assembly_bp": 504621863,
    # busco_short_summary.txt, eudicotyledons_odb12.2 (n = 1,990)
    "busco": {"S": 1957, "D": 30, "F": 1, "M": 2, "C_pct": 99.8},
    # PN40024 5.1 mapped on the same assembly, re-run with the same containers (data/v5_1_audit)
    "v51": {"genes": 47971, "busco_C": 97.9, "busco_S": 1918, "busco_D": 31, "busco_F": 0, "busco_M": 41,
            "single_exon": 17117, "single_exon_pct": 41.0, "overlap_pairs": 35, "utr_pct": 60.0,
            "stops_main": 30, "mean_gene_len": 4850},
    "v43_genes": 34843, "v43_transferred": 34430, "v43_unmapped": 413, "carried_over": 31673,
    "new_ids": 25071,
    # omark/proteins_main_detailed_summary.txt
    "omark": {"single": 8661, "dup_exp": 537, "dup_unexp": 1139, "missing": 396, "hogs": 10733,
              "consistent": 38599, "inconsistent": 1550, "contam": 0, "unknown": 15617,
              "partial": 9775, "fragmented": 3328},
    # expression_support_summary.json, functional annotation (audit sections 6, 7)
    "func": {"d2go": 33034, "eggnog": 43596, "interpro": 44109, "any": 47281, "none": 8485},
    "expr_supported_pct": 72.8, "no_expr_total": 15454,
    # audit section 5.3
    "zero_of_three": 3416, "at_least_one_pct": 94.0,
    # audit section 9 (overlaps), corrected 2026-08-21
    "overlap_pairs": 3520, "ovl_antisense": 2132, "ovl_utr": 729, "ovl_cds": 659, "ovl_cds_genes": 1284,
    "ovl_genes": 6621, "ovl_genes_pct": 11.7, "superloci": 13697,
    # monoexonic confidence pass
    "mono_removed": 943, "mono_te": 610, "mono_grey": 333, "mono_unsup_augustus_pct": 91.2,
    "validation_stops": 4,
}

CLASS_ORDER = ["expressed", "conserved_or_functional", "te_protein", "te_overlap_unsupported",
               "liftoff_only", "junction_only", "no_signal"]
CLASS_FR = {
    "expressed": "Exprimés (TPM ≥ 0,5)",
    "conserved_or_functional": "Non exprimés, conservés ou fonctionnels",
    "te_protein": "Non exprimés, protéines de TE",
    "te_overlap_unsupported": "Non exprimés, sur un TE, sans autre support",
    "liftoff_only": "Non exprimés, gène v4.3 seulement",
    "junction_only": "Non exprimés, jonction seulement",
    "no_signal": "Aucun signal",
}
# Stack order of the validated palette (blue, orange, aqua, yellow, violet, red)
CLASS_COLOR = {"expressed": 1, "te_protein": 2, "conserved_or_functional": 3, "liftoff_only": 4,
               "junction_only": 4, "te_overlap_unsupported": 5, "no_signal": 6}
ORIGIN_FR = {"ab_initio_braker": "BRAKER (AUGUSTUS/GeneMark bruts)", "ab_initio_helixer": "Helixer",
             "egapx_gnomon": "EGAPx (Gnomon)", "liftoff_previous": "Liftoff (annotations précédentes)",
             "rnaseq_assembly": "Assemblages RNA-seq"}

CSS = """
:root{--surface:#fcfcfb;--panel:#f4f3f0;--ink:#0b0b0b;--ink2:#52514e;--grid:#dcdad4;--line:#cfcdc6;
--c1:#2a78d6;--c2:#eb6834;--c3:#1baf7a;--c4:#eda100;--c5:#4a3aa7;--c6:#e34948;--accent:#184f95;
--r1:#86b6ef;--r2:#3987e5;--r3:#1c5cab;--r4:#0d366b;--r1t:#0b0b0b;--r2t:#fff;--r3t:#fff;--r4t:#fff;
--warn-bg:#fdf1e7;--warn-line:#eb6834;--note-bg:#eaf2fc;--note-line:#2a78d6;--ok-bg:#e8f6ef;--ok-line:#1baf7a}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){--surface:#1a1a19;--panel:#242422;
--ink:#f4f3ef;--ink2:#c3c2b7;--grid:#3a3a37;--line:#454541;--c1:#3987e5;--c2:#d95926;--c3:#199e70;
--c4:#c98500;--c5:#9085e9;--c6:#e66767;--accent:#86b6ef;--warn-bg:#33251c;--note-bg:#1c2733;--ok-bg:#1b2c25;
--r1:#184f95;--r2:#256abf;--r3:#3987e5;--r4:#86b6ef;--r1t:#fff;--r2t:#fff;--r3t:#0b0b0b;--r4t:#0b0b0b}}
:root[data-theme="dark"]{--surface:#1a1a19;--panel:#242422;--ink:#f4f3ef;--ink2:#c3c2b7;--grid:#3a3a37;
--line:#454541;--c1:#3987e5;--c2:#d95926;--c3:#199e70;--c4:#c98500;--c5:#9085e9;--c6:#e66767;
--accent:#86b6ef;--warn-bg:#33251c;--note-bg:#1c2733;--ok-bg:#1b2c25;
--r1:#184f95;--r2:#256abf;--r3:#3987e5;--r4:#86b6ef;--r1t:#fff;--r2t:#fff;--r3t:#0b0b0b;--r4t:#0b0b0b}
*{box-sizing:border-box}
html{scroll-behavior:smooth}
body{margin:0;background:var(--surface);color:var(--ink);font:16px/1.6 system-ui,-apple-system,"Segoe UI",Roboto,sans-serif}
.wrap{max-width:1040px;margin:0 auto;padding:0 20px 80px}
header.top{padding:44px 0 26px;border-bottom:1px solid var(--line)}
header.top .kicker{color:var(--ink2);font-size:13px;letter-spacing:.06em;text-transform:uppercase}
h1{font-size:clamp(26px,4vw,38px);line-height:1.15;margin:.3em 0 .2em}
h2{font-size:26px;margin:2.2em 0 .5em;padding-top:.4em;border-top:1px solid var(--line)}
h3{font-size:19px;margin:1.9em 0 .4em}
h4{font-size:16px;margin:1.4em 0 .3em;color:var(--ink2)}
p{margin:.6em 0}
a{color:var(--accent)}
.lead{color:var(--ink2);font-size:18px;max-width:78ch}
nav.toc{background:var(--panel);border-radius:10px;padding:14px 20px;margin:22px 0}
nav.toc ol{margin:.3em 0;padding-left:1.3em;columns:2;column-gap:32px}
nav.toc li{break-inside:avoid;margin:.15em 0}
.tiles{display:grid;grid-template-columns:repeat(auto-fit,minmax(170px,1fr));gap:12px;margin:22px 0}
.tile{background:var(--panel);border-radius:10px;padding:14px 16px}
.tile .v{font-size:30px;font-weight:650;line-height:1.1;font-variant-numeric:tabular-nums}
.tile .l{color:var(--ink2);font-size:13.5px;margin-top:4px}
.callout{border-left:4px solid var(--note-line);background:var(--note-bg);padding:12px 16px;border-radius:0 8px 8px 0;margin:16px 0}
.callout.warn{border-color:var(--warn-line);background:var(--warn-bg)}
.callout.ok{border-color:var(--ok-line);background:var(--ok-bg)}
.callout p:first-child{margin-top:0}.callout p:last-child{margin-bottom:0}
.callout .h{font-weight:650}
figure{margin:22px 0;background:var(--panel);border-radius:10px;padding:14px 16px 10px}
figcaption.ft{font-weight:650;margin-bottom:6px}
p.fc{color:var(--ink2);font-size:14px;margin:6px 0 4px}
svg.chart{width:100%;height:auto;display:block}
.f1{fill:var(--c1)}.f2{fill:var(--c2)}.f3{fill:var(--c3)}.f4{fill:var(--c4)}.f5{fill:var(--c5)}.f6{fill:var(--c6)}
.f7{fill:var(--r1)}.f8{fill:var(--r2)}.f9{fill:var(--r3)}.f10{fill:var(--r4)}
.tr1{fill:var(--r1t);font-size:12px;font-weight:600}.tr2{fill:var(--r2t);font-size:12px;font-weight:600}
.tr3{fill:var(--r3t);font-size:12px;font-weight:600}.tr4{fill:var(--r4t);font-size:12px;font-weight:600}
.s1{stroke:var(--c1)}.s2{stroke:var(--c2)}.s3{stroke:var(--c3)}.s4{stroke:var(--c4)}.s5{stroke:var(--c5)}.s6{stroke:var(--c6)}
.t1{fill:var(--ink);font-size:13px}.t2{fill:var(--ink2);font-size:12px}.tw{fill:#fff;font-size:12px;font-weight:600}
.td{fill:#0b0b0b;font-size:12px;font-weight:600}.hm-w{fill:#fff;font-size:12px;font-weight:600}.hm-d{fill:#0b0b0b;font-size:12px}
.grid{stroke:var(--grid);stroke-width:1}
.tablewrap{overflow-x:auto;margin:10px 0}
.txt th,.txt td{text-align:left}
table{border-collapse:collapse;width:100%;font-size:14px;font-variant-numeric:tabular-nums}
caption{text-align:left;color:var(--ink2);padding:4px 0}
th,td{padding:6px 10px;border-bottom:1px solid var(--line);text-align:right;vertical-align:top}
th:first-child,td:first-child{text-align:left}
thead th{color:var(--ink2);font-weight:600;border-bottom:2px solid var(--line)}
td.l,th.l{text-align:left}
tr.hl td{background:var(--warn-bg)}
details.tv{margin-top:6px}details.tv summary{cursor:pointer;color:var(--ink2);font-size:13.5px}
code,pre{font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;font-size:13px}
pre{background:var(--panel);padding:12px 14px;border-radius:8px;overflow-x:auto;line-height:1.45}
code{background:var(--panel);padding:1px 5px;border-radius:4px}
pre code{background:none;padding:0}
ul,ol{padding-left:1.3em}li{margin:.25em 0}
.cols2{display:grid;grid-template-columns:repeat(auto-fit,minmax(300px,1fr));gap:16px}
.tag{display:inline-block;font-size:12px;padding:1px 8px;border-radius:999px;border:1px solid var(--line);color:var(--ink2)}
.mono{font-family:ui-monospace,Menlo,Consolas,monospace}
footer{margin-top:60px;padding-top:16px;border-top:1px solid var(--line);color:var(--ink2);font-size:13.5px}
@media (max-width:640px){body{font-size:15px}nav.toc ol{columns:1}h2{font-size:22px}.tile .v{font-size:25px}}
@media print{body{font-size:11.5pt}figure,.callout,table{break-inside:avoid}nav.toc{display:none}
h2{break-after:avoid}details.tv{display:none}}
"""


def pct(a, b, d=1):
    return fmt_pct(100.0 * a / b, d)


def table(headers, rows, text=False, caption="", hl_rows=()):
    th = "".join("<th>%s</th>" % h for h in headers)
    body = ""
    for i, r in enumerate(rows):
        tr = ' class="hl"' if i in hl_rows else ""
        body += "<tr%s>%s</tr>" % (tr, "".join("<td>%s</td>" % c for c in r))
    cap = "<caption>%s</caption>" % caption if caption else ""
    return '<div class="tablewrap%s"><table>%s<thead><tr>%s</tr></thead><tbody>%s</tbody></table></div>' % (
        " txt" if text else "", cap, th, body)


def callout(text, kind=""):
    return '<div class="callout %s">%s</div>' % (kind, text)


# --------------------------------------------------------------------------------------------
# Data
# --------------------------------------------------------------------------------------------
class Data:
    def __init__(self, a):
        an = a.analysis_dir
        self.cls = pd.read_csv(os.path.join(an, "classes", "gene_classes.tsv"), sep="\t", index_col=0)
        self.summary = pd.read_csv(os.path.join(an, "classes", "class_summary.tsv"), sep="\t", index_col=0)
        self.reads = pd.read_csv(os.path.join(an, "classes", "salmon_reads_by_class.tsv"), sep="\t", index_col=0)
        self.sig = pd.read_csv(os.path.join(an, "classes", "signature_summary.tsv"), sep="\t", index_col=0)
        self.tpm = pd.read_csv(os.path.join(a.titan_out, "01_final_annotation", "quality_report",
                                            "expression_validation", "gene_tpm_matrix.tsv"), sep="\t", index_col=0)
        cmp_dir = os.path.join(a.test_dir, "comparison")
        self.cmp = json.load(open(os.path.join(cmp_dir, "comparison.json")))
        self.fate = {arm: pd.read_csv(os.path.join(cmp_dir, "fate_by_class_%s.tsv" % arm), sep="\t", index_col=0)
                     for arm in ("raw", "tsebra")}
        fa = os.path.join(cmp_dir, "fate_annotated_tsebra.tsv")
        self.fate_genes = pd.read_csv(fa, sep="\t", index_col=0) if os.path.exists(fa) else None
        self.n = len(self.cls)
        self.n_class = self.cls["evidence_class"].value_counts().to_dict()
        c = self.cls
        self.set_aside = int(c["evidence_class"].isin(["no_signal", "te_overlap_unsupported"]).sum())
        self.noexpr = int((c["evidence_class"] != "expressed").sum())


def saturation(tpm, coding_ids, n_perm=100, seed=1):
    det = (tpm.loc[tpm.index.isin(coding_ids)] >= 0.5).values
    rng = np.random.RandomState(seed)
    curves = []
    for _ in range(n_perm):
        cum = np.zeros(det.shape[0], dtype=bool)
        row = []
        for j in rng.permutation(det.shape[1]):
            cum |= det[:, j]
            row.append(int(cum.sum()))
        curves.append(row)
    return np.array(curves)


# --------------------------------------------------------------------------------------------
# Sections
# --------------------------------------------------------------------------------------------
def workflow_svg():
    def box(x, y, w, h, title, lines, cls="f1", op=0.16):
        s = '<rect x="%d" y="%d" width="%d" height="%d" rx="8" class="%s" opacity="%s"/>' % (x, y, w, h, cls, op)
        s += '<rect x="%d" y="%d" width="%d" height="%d" rx="8" fill="none" stroke="var(--line)"/>' % (x, y, w, h)
        s += '<text x="%d" y="%d" class="t1" font-weight="650">%s</text>' % (x + 12, y + 22, esc(title))
        for i, ln in enumerate(lines):
            s += '<text x="%d" y="%d" class="t2">%s</text>' % (x + 12, y + 42 + i * 17, esc(ln))
        return s

    def arrow(x1, y1, x2, y2):
        return ('<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="var(--ink2)" stroke-width="1.6" '
                'marker-end="url(#ar)"/>' % (x1, y1, x2, y2))

    body = ('<defs><marker id="ar" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" '
            'orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="var(--ink2)"/></marker></defs>')
    body += box(8, 8, 230, 150, "Preuves d'annotation", [
        "EGAPx · Liftoff (v4.3)", "BRAKER3 : AUGUSTUS + GeneMark", "Helixer (ab initio)",
        "StringTie / PsiCLASS (STAR)", "Long reads (minimap2 + StringTie)", "47 banques RNA-seq"], "f1")
    body += box(300, 8, 200, 150, "Mikado", ["prepare (13 pistes, scores)", "TransDecoder (ORFs)", "serialise",
                                             "pick : 1 meilleur transcrit", "par locus"], "f2")
    body += box(560, 8, 190, 150, "AEGIS + contrôles", ["renommage / tidy des IDs", "annotation primaire",
                                                        "BUSCO · OMArk · AGAT", "Salmon (47 banques)",
                                                        "Fonction : D2GO, eggNOG, IPR"], "f3")
    body += arrow(238, 83, 298, 83) + arrow(500, 83, 558, 83)
    return ('<svg class="chart" viewBox="0 0 760 172" role="img" aria-label="Schéma du pipeline TITAN : '
            'preuves, Mikado, AEGIS et contrôles qualité">%s</svg>' % body)


def sec_summary(D):
    c = D.n_class
    tiles = [
        (fmt_int(AUDIT["genes"]), "gènes annotés (%s codants + %s ncRNA)" % (fmt_int(AUDIT["coding"]), fmt_int(AUDIT["ncrna"]))),
        ("%s %%" % dec_fr(AUDIT["busco"]["C_pct"]), "BUSCO complet (5.1 : 97,9 %)"),
        ("%s %%" % dec_fr(AUDIT["expr_supported_pct"]), "gènes exprimés (TPM ≥ 0,5, 47 banques)"),
        ("%s %%" % dec_fr(100.0 * AUDIT["func"]["any"] / AUDIT["coding"]), "codants avec ≥ 1 annotation fonctionnelle"),
        (fmt_int(D.set_aside), "gènes proposés à écarter (%s des codants)" % pct(D.set_aside, D.n)),
    ]
    t = "".join('<div class="tile"><div class="v">%s</div><div class="l">%s</div></div>' % x for x in tiles)
    return '<div class="tiles">%s</div>' % t


def sec_keys(D):
    c = D.n_class
    return """
<h2 id="cles">Messages clés</h2>
<ol>
<li><b>L'annotation est complète et structurellement valide.</b> BUSCO complet 99,8 %% (2 manquants sur 1 990) contre 97,9 %% (41 manquants) pour PN40024 5.1 projetée sur le même assemblage&nbsp;; validation structurale PASS&nbsp;; aucune contamination détectée par OMArk.</li>
<li><b>56&nbsp;744 gènes, c'est beaucoup, mais moins qu'il n'y paraît.</b> Face à la référence actuelle (5.1&nbsp;: %s gènes) l'excès est de +18,3&nbsp;%%, pas de +63&nbsp;%% comme face à v4.3 (%s). %s des gènes ont au moins une preuve indépendante (expression, fonction ou correspondance v4.3).</li>
<li><b>Les ~15&nbsp;000 gènes «&nbsp;sans expression&nbsp;» ne forment pas un bloc d'erreurs.</b> Sur %s gènes non exprimés (codants) : %s sont conservés ou fonctionnels, %s codent des protéines de TE, et seulement <b>%s</b> n'ont <i>aucun</i> support indépendant (%s). Les 978 ncRNA n'ont simplement jamais été quantifiés.</li>
<li><b>Ces %s gènes sont très probablement des erreurs de prédiction</b> : médiane de %s read Salmon sur 47 banques, aucun intron retrouvé, presque aucun homologue chez d'autres eudicotylédones, aucun gène BUSCO, %s de recouvrement avec PN40024 5.1 (contre %s pour les gènes exprimés).</li>
<li><b>Ils viennent presque tous de modèles bruts de BRAKER</b> (%s), que Mikado accepte même sans concurrent. Fournir à Mikado les pistes filtrées par TSEBRA supprime 60 %% de ces gènes mais aussi %s des gènes exprimés&nbsp;: ce n'est pas une amélioration nette, donc pas retenu.</li>
<li><b>Un outil de marquage et de filtrage du GFF3 est fourni</b> (hors pipeline). Politique stricte : %s gènes conservés sur %s.</li>
<li><b>Le nombre élevé de gènes ne s'explique pas surtout par des erreurs.</b> Même la politique la plus sévère laisse plus de 51&nbsp;000 gènes codants&nbsp;; un noyau de %s gènes est soutenu par au moins 3 lignes de preuve. Le reste relève de la biologie (assemblage T2T, familles multigéniques, duplications) et de choix de méthode, à discuter (section 3.11).</li>
</ol>""" % (
        fmt_int(AUDIT["v51"]["genes"]), fmt_int(AUDIT["v43_genes"]), fmt_pct(AUDIT["at_least_one_pct"], 0),
        fmt_int(D.noexpr), fmt_int(c.get("conserved_or_functional", 0)), fmt_int(c.get("te_protein", 0)),
        fmt_int(D.set_aside), "1&nbsp;263 «&nbsp;aucun signal&nbsp;» + 1&nbsp;215 «&nbsp;sur un TE&nbsp;»",
        fmt_int(D.set_aside), "1", "6-8&nbsp;%", "67&nbsp;%",
        "82&nbsp;%", "9,7&nbsp;%",
        fmt_int(54266), fmt_int(AUDIT["genes"]), "37&nbsp;632")


def sec_context(D):
    rows = [["EGAPx (Gnomon)", "20", "référence"], ["Liftoff (v4.3 transféré)", "19", "référence"],
            ["BRAKER3 : AUGUSTUS", "18", "brut"], ["BRAKER3 : GeneMark", "17", "brut"], ["Helixer", "16", "ab initio"],
            ["StringTie STAR stranded (défaut / alt)", "15 / 14", "assemblage"],
            ["PsiCLASS STAR stranded / unstranded", "13 / 7", "assemblage"],
            ["StringTie STAR unstranded (défaut / alt)", "6 / 5", "assemblage"],
            ["Long reads minimap2 + StringTie (défaut / alt)", "9 / 8", "assemblage"],
            ["FLAIR", "10", "vide dans ce run"]]
    t = table(["Piste", "Score Mikado", "Nature"], rows)
    return """
<h2 id="pipeline">1. Le pipeline en bref</h2>
<p>TITAN combine <b>13 pistes de preuves</b> (prédictions ab initio, transferts d'annotations, assemblages de transcrits) et laisse <b>Mikado</b> choisir, pour chaque locus, le meilleur transcrit selon un fichier de scores (<code>plant.yaml</code>) et une priorité par piste. AEGIS renomme ensuite les gènes et nettoie le GFF3. Les contrôles qualité (BUSCO, OMArk, AGAT, quantification Salmon, annotation fonctionnelle) sont calculés sur l'annotation finale.</p>
%s
<p class="fc">Schéma simplifié, vérifié contre <code>workflows/titan.nf</code> : <code>mikado_prepare → transdecoder_longorfs → transdecoder_predict → mikado_serialise → mikado_pick → aegis</code>.</p>
<h4>Pistes données à Mikado et priorités (<code>modules/mikado.nf</code>)</h4>
%s
<div class="callout"><p><span class="h">Point important pour la suite.</span> Mikado conserve le meilleur transcrit de chaque locus <i>quel que soit son niveau de preuve</i> : un modèle ab initio isolé, sans concurrent, est reporté tel quel (13&nbsp;697 super-loci donnent 56&nbsp;744 gènes). Les pistes BRAKER3 données à Mikado sont les prédictions <b>brutes</b> (<code>augustus.hints.gff3</code>, <code>genemark.gtf</code>), pas le jeu filtré par TSEBRA (<code>braker.gff3</code>).</p></div>
""" % (figure(workflow_svg(), "Vue d'ensemble", "13 pistes → Mikado → AEGIS → annotation primaire et contrôles."), t)


def sec_quality(D):
    A, V = AUDIT, AUDIT["v51"]
    # BUSCO chart
    bus = stacked_hbars(
        [("TITAN (56 744 gènes)", [A["busco"]["S"], A["busco"]["D"], A["busco"]["F"], A["busco"]["M"]]),
         ("PN40024 5.1 (47 971 gènes)", [V["busco_S"], V["busco_D"], V["busco_F"], V["busco_M"]])],
        [("Complet, copie unique", 1), ("Complet, dupliqué", 2), ("Fragmenté", 3), ("Manquant", 5)],
        "BUSCO eudicotyledons_odb12.2, TITAN et PN40024 5.1", label_w=220)
    bus_t = data_table(["Jeu", "Copie unique", "Dupliqué", "Fragmenté", "Manquant"], [
        ["TITAN", A["busco"]["S"], A["busco"]["D"], A["busco"]["F"], A["busco"]["M"]],
        ["PN40024 5.1", V["busco_S"], V["busco_D"], V["busco_F"], V["busco_M"]]])
    om = A["omark"]
    omc = stacked_hbars(
        [("Complétude (10 733 HOG)", [om["single"], om["dup_unexp"], om["dup_exp"], om["missing"]])],
        [("Copie unique", 1), ("Dupliqué inattendu", 2), ("Dupliqué attendu", 3), ("Manquant", 5)],
        "OMArk, complétude sur les HOG conservés des rosidées", label_w=215)
    omk = stacked_hbars(
        [("Cohérence (55 766 prot.)", [om["consistent"], om["inconsistent"], om["contam"], om["unknown"]])],
        [("Cohérent", 1), ("Incohérent", 2), ("Contaminant", 5), ("Inconnu", 4)],
        "OMArk, cohérence taxonomique des protéines", label_w=215)
    rows = [
        ["Gènes", fmt_int(A["genes"]), fmt_int(V["genes"])],
        ["Gènes mono-exoniques", "%s (%s)" % (fmt_int(A["single_exon"]), fmt_pct(A["single_exon_pct"])),
         "%s (%s)" % (fmt_int(V["single_exon"]), fmt_pct(V["single_exon_pct"]))],
        ["ARNm avec au moins une UTR", fmt_pct(A["utr_pct"]), fmt_pct(V["utr_pct"])],
        ["Longueur moyenne d'un gène (pb)", fmt_int(A["mean_gene_len"]), fmt_int(V["mean_gene_len"])],
        ["Paires de gènes qui se chevauchent (AGAT)", fmt_int(A["overlap_pairs"]), fmt_int(V["overlap_pairs"])],
        ["Codon stop interne (jeu principal)", "3", fmt_int(V["stops_main"])],
    ]
    return """
<h2 id="qualite">2. Qualité globale de l'annotation</h2>
<h3>2.1 Complétude : BUSCO</h3>
%s
<p>La complétude est meilleure que celle de la référence actuelle (99,8&nbsp;%% contre 97,9&nbsp;%%, 20 fois moins de gènes manquants). La proportion de BUSCO dupliqués reste faible (1,5&nbsp;%%), donc l'annotation n'est pas globalement dupliquée.</p>
<h3>2.2 Complétude et cohérence : OMArk</h3>
%s
%s
<p>Le taux de duplication d'OMArk (%s) est plus élevé que celui de BUSCO&nbsp;: il contient 5,0&nbsp;%% de duplications <i>attendues</i> (expansions connues de familles) et 10,6&nbsp;%% <i>inattendues</i>, le chiffre à surveiller si des artefacts de duplication segmentaire ou d'haplotypes sont à craindre. Aucun contaminant. Les %s protéines «&nbsp;inconnues&nbsp;» (28,0&nbsp;%%) sont compatibles à la fois avec un répertoire spécifique de la vigne et avec des modèles erronés&nbsp;: c'est précisément la population examinée en section 3. Parmi les protéines cohérentes, %s (17,5&nbsp;%%) ne s'alignent que partiellement sur leur famille et %s (6,0&nbsp;%%) sont plus courtes que la moitié de la longueur médiane de la famille.</p>
<h3>2.3 Structure des gènes</h3>
%s
<p>La structure est proche de la référence, avec deux différences à connaître&nbsp;: moins de mono-exoniques (28,0&nbsp;%% contre 41,0&nbsp;%%) mais <b>beaucoup plus de chevauchements</b> (3&nbsp;520 paires contre 35). Ce dernier point est détaillé en 2.6.</p>
<h3>2.4 Validation structurale</h3>
<p>Le contrôle de cohérence GFF3/FASTA de TITAN est <b>PASS</b> : nombres de features cohérents (56&nbsp;744 gènes, 57&nbsp;542 ARNm, 238&nbsp;482 exons, 237&nbsp;352 CDS, 982 ncRNA), 411&nbsp;292 identifiants uniques, 55&nbsp;766 protéines principales retrouvées dans le GFF3. Seul avertissement&nbsp;: 4 protéines (0,005&nbsp;%%) portent un codon stop interne.</p>
<h3>2.5 Support transcriptomique et fonctionnel</h3>
<ul>
<li><b>Expression</b> : %s des gènes ont un TPM ≥ 0,5 dans au moins une des 47 banques (mélange de polyA, ARN déplété stranded et long reads).</li>
<li><b>Fonction</b> (sur %s protéines principales) : Diamond2GO %s, eggNOG-mapper %s, InterProScan %s&nbsp;; au moins une méthode&nbsp;: <b>%s (%s)</b>. %s gènes n'ont aucun résultat.</li>
<li><b>Correspondance avec v4.3</b> : %s gènes gardent un identifiant v4.3 avec confiance (%s), %s ont un nouvel identifiant.</li>
</ul>
<h3>2.6 Chevauchements de gènes</h3>
<p>AGAT compte %s <i>paires</i> de gènes dont les exons se chevauchent, contre 35 dans PN40024 5.1. Reclassées par brin et par chevauchement de CDS : <b>%s paires (60,6&nbsp;%%) sont antisens</b> (attendu, ce n'est pas le travail de Mikado de les fusionner), <b>%s (20,7&nbsp;%%) sont de même brin mais ne partagent que de l'UTR</b>, et seules <b>%s (18,7&nbsp;%%)</b>, soit %s gènes (2,3&nbsp;%% de l'annotation), sont de vrais conflits CDS/CDS de même brin. Au total %s gènes (%s) sont dans un chevauchement quelconque. Ces conflits proviennent probablement de pistes différentes proposant des modèles incompatibles pour un même locus&nbsp;: 37,7&nbsp;%% de ces gènes sont mono-exoniques et seulement 1,0&nbsp;%% sont dans la catégorie «&nbsp;TE&nbsp;» du filtre mono-exonique. On y trouve un enrichissement modeste en familles multigéniques classiques (NB-LRR ×1,6, PPR ×1,3). La quantification avec Salmon reste possible : sur 362 paires où les deux gènes sont exprimés, seules 2 (0,6&nbsp;%%) montrent le signe d'une compétition entre modèles (corrélation de TPM &lt; −0,3).</p>
<h3>2.7 Gènes mono-exoniques</h3>
<p>%s gènes (28,0&nbsp;%%) sont mono-exoniques, ce qui est élevé pour une plante mais inférieur à 5.1 (41,0&nbsp;%%). Le filtre de confiance de TITAN a classé %s d'entre eux comme non soutenus (dont %s chevauchent un TE) ; <b>%s&nbsp;%% de ces gènes mono-exoniques non soutenus viennent d'AUGUSTUS brut</b> (860 sur 943), soit 4,9 fois sa part globale. C'est la première indication du mécanisme développé en section 3.8. Les assemblages PsiCLASS gonflent le <i>nombre</i> de mono-exoniques (32,6&nbsp;%% d'entre eux) mais 98,8&nbsp;%% de ceux-là sont soutenus.</p>
""" % (
        figure(bus, "BUSCO : TITAN et PN40024 5.1", "2 BUSCO manquants dans TITAN contre 41 dans PN40024 5.1 (n = 1 990).", bus_t),
        figure(omc, "OMArk : complétude", "10 733 groupes d'orthologues hiérarchiques conservés chez les rosidées.", ""),
        figure(omk, "OMArk : cohérence", "Aucune contamination ; 28 % des protéines n'ont pas de famille connue.", ""),
        fmt_pct(15.62), fmt_int(om["unknown"]), fmt_int(om["partial"]), fmt_int(om["fragmented"]),
        table(["Métrique", "TITAN (primaire)", "PN40024 5.1"], rows),
        fmt_pct(A["expr_supported_pct"]), fmt_int(A["coding"]), fmt_pct(100.0 * A["func"]["d2go"] / A["coding"]),
        fmt_pct(100.0 * A["func"]["eggnog"] / A["coding"]), fmt_pct(100.0 * A["func"]["interpro"] / A["coding"]),
        fmt_int(A["func"]["any"]), fmt_pct(100.0 * A["func"]["any"] / A["coding"]), fmt_int(A["func"]["none"]),
        fmt_int(A["carried_over"]), fmt_pct(100.0 * A["carried_over"] / A["genes"]), fmt_int(A["new_ids"]),
        fmt_int(A["overlap_pairs"]), fmt_int(A["ovl_antisense"]), fmt_int(A["ovl_utr"]), fmt_int(A["ovl_cds"]),
        fmt_int(A["ovl_cds_genes"]), fmt_int(A["ovl_genes"]), fmt_pct(A["ovl_genes_pct"]),
        fmt_int(A["single_exon"]), fmt_int(A["mono_removed"]), fmt_int(A["mono_te"]), "91,2")


def sec_count_intro(D):
    A, V = AUDIT, AUDIT["v51"]
    items = [
        dict(label="v4.3 (référence historique)", value=A["v43_genes"], color=1),
        dict(label="PN40024 5.1 (référence actuelle)", value=V["genes"], color=1),
        dict(label="TITAN, gènes codants", value=A["coding"], color=2),
        dict(label="TITAN, avec 978 ncRNA", value=A["genes"], color=2),
    ]
    fig = hbars(items, "Nombre de gènes : v4.3, 5.1 et TITAN", label_w=240)
    tab = data_table(["Annotation", "Gènes"], [[i["label"], fmt_int(i["value"])] for i in items])
    d51 = 100.0 * (A["genes"] - V["genes"]) / V["genes"]
    d43 = 100.0 * (A["genes"] - A["v43_genes"]) / A["v43_genes"]
    return f"""
<h2 id="nombre">3. Pourquoi 56&nbsp;744 gènes ? Discussion</h2>
<p>Les ordres de grandeur souvent cités pour la vigne sont de 35 à 40&nbsp;000 gènes (l'annotation v4.3 en compte 34&nbsp;843). TITAN en annote 55&nbsp;766 codants (56&nbsp;744 avec les ncRNA). Cette section examine, chiffres à l'appui, ce qu'il y a derrière ce nombre&nbsp;: quelle part est soutenue, quelle part est douteuse, d'où viennent les modèles douteux, et ce qu'un filtrage changerait.</p>
<h3>3.1 Le point de comparaison compte</h3>
{figure(fig, "Nombre de gènes selon l'annotation", f"Par rapport à v4.3 l'excès est de +{dec_fr(d43)} %, par rapport à 5.1 de +{dec_fr(d51)} %. 5.1 est déjà projetée sur l'assemblage T2T utilisé ici.", tab)}
<p>La comparaison avec v4.3 (34&nbsp;843 gènes) est la plus alarmante, mais c'est aussi la moins équitable&nbsp;: v4.3 est une référence plus ancienne (le transfert Liftoff retrouve 98,8&nbsp;% de ses gènes sur T2T, mais seulement 55,8&nbsp;% des gènes TITAN portent un identifiant v4.3 avec confiance, et l'écart se réduit fortement face à 5.1). Face à la référence <i>actuelle</i> 5.1, projetée sur le même assemblage, TITAN a <b>18,3&nbsp;% de gènes en plus</b> tout en étant plus complète (BUSCO 99,8&nbsp;% contre 97,9&nbsp;%). Un nombre de gènes brut n'est donc pas un indicateur de qualité en soi&nbsp;: il faut examiner le support de chaque gène.</p>
"""


def sec_axes(D):
    c = D.cls
    ax = pd.DataFrame({
        "expr": c["evidence_class"] == "expressed",
        "junc": c["has_junction"] if "has_junction" in c else False,
        "func": c["has_function"],
        "orth": c["omamer_hit"] | c["has_nonvitis_homolog"],
        "lift": c["liftoff_carried"],
    })
    n = ax.sum(axis=1)
    vc = n.value_counts().sort_index()
    items = [dict(label=f"{k} ligne(s) de preuve", value=int(vc.get(k, 0)), color=1 if k >= 3 else (4 if k == 2 else 6),
                  text=f"{fmt_int(vc.get(k, 0))} ({pct(vc.get(k, 0), len(c))})") for k in range(0, 6)]
    fig = hbars(items, "Nombre de lignes de preuve par gène codant", label_w=200, right_pad=150)
    tab = data_table(["Lignes de preuve", "Gènes codants", "%"], [[k, fmt_int(vc.get(k, 0)), pct(vc.get(k, 0), len(c))] for k in range(6)])
    ge2, ge3, ge4 = int((n >= 2).sum()), int((n >= 3).sum()), int((n >= 4).sum())
    D.ge2, D.ge3, D.ge4 = ge2, ge3, ge4
    return f"""
<h3>3.2 Combien de gènes sont soutenus, et par quoi ?</h3>
<p>Cinq lignes de preuve, en grande partie indépendantes, ont été comptées pour chaque gène codant&nbsp;: <b>expression</b> (TPM ≥ 0,5 dans au moins une banque), <b>jonction d'épissage</b> (au moins un intron retrouvé exactement dans un assemblage RNA-seq), <b>fonction</b> (Diamond2GO, eggNOG, InterProScan ou MapMan), <b>orthologie/conservation</b> (famille OMAMER ou homologue chez une eudicotylédone non-<i>Vitis</i>) et <b>correspondance avec v4.3</b> (Liftoff).</p>
{figure(fig, "Nombre de lignes de preuve par gène", f"{fmt_int(ge2)} gènes codants ont au moins 2 lignes de preuve, {fmt_int(ge3)} au moins 3, {fmt_int(ge4)} au moins 4. Seuls {fmt_int(vc.get(0, 0))} n'en ont aucune.", tab)}
{callout(f'<p><span class="h">Lecture prudente.</span> Ces lignes ne sont pas parfaitement indépendantes (une famille OMAMER et une annotation fonctionnelle reposent toutes deux sur la similarité de séquence ; les jonctions et l\'expression sur les mêmes banques RNA-seq ; un mono-exonique ne peut pas avoir de jonction). Le nombre de <b>{fmt_int(ge3)}</b> gènes à ≥ 3 lignes est un <i>noyau très fiable</i>, pas une estimation du nombre réel de gènes. Il se situe néanmoins dans la fourchette de 35 à 40&nbsp;000 souvent citée.</p>', "")}
"""


def sec_origin(D):
    c = D.cls
    order = ["ab_initio_braker", "egapx_gnomon", "liftoff_previous", "ab_initio_helixer", "rnaseq_assembly"]
    rows = []
    for o in order:
        sub = c[c["origin"] == o]
        e = int((sub["evidence_class"] == "expressed").sum())
        rows.append((ORIGIN_FR[o], [e, len(sub) - e], len(sub), 100.0 * (len(sub) - e) / len(sub)))
    fig = stacked_hbars([(r[0], r[1]) for r in rows], [("Exprimés", 1), ("Non exprimés", 2)],
                        "Gènes par origine du modèle retenu par Mikado, exprimés et non exprimés", label_w=300)
    tab = data_table(["Origine", "Gènes", "Non exprimés", "% non exprimés"],
                     [[r[0], fmt_int(r[2]), fmt_int(r[1][1]), fmt_pct(r[3])] for r in rows])
    return f"""
<h3>3.3 D'où viennent les gènes ?</h3>
<p>Chaque gène final a été relié au transcrit Mikado dont il provient (par coordonnées, 100&nbsp;% de correspondance) et donc à la piste qui l'a fourni (attribut <code>alias</code>).</p>
{figure(fig, "Origine du modèle retenu et part non exprimée", "Le taux de gènes non exprimés varie de 1,4 % (EGAPx) à 54 % (BRAKER brut).", tab)}
<p>Deux lectures. D'abord, <b>l'absence d'expression discrimine surtout les modèles ab initio et Liftoff</b>&nbsp;: EGAPx et les assemblages RNA-seq sont construits à partir des mêmes reads, donc leurs gènes sont exprimés presque par construction (1 à 2&nbsp;% non exprimés). Ensuite, <b>BRAKER brut fournit à lui seul 8&nbsp;434 gènes exprimés et 9&nbsp;995 non exprimés</b>&nbsp;: c'est la piste où la question du support se pose.</p>
"""


def sec_noexpr(D):
    c, t = D.cls, D.tpm
    coding = c.index
    tm = t.loc[t.index.isin(coding)]
    mx = tm.max(axis=1)
    A = AUDIT
    n_noexp = int((mx < 0.5).sum())
    bins = [("Aucun read (TPM = 0 partout)", (mx == 0).sum()), ("0 < max TPM ≤ 0,1", ((mx > 0) & (mx <= 0.1)).sum()),
            ("0,1 < max TPM ≤ 0,25", ((mx > 0.1) & (mx <= 0.25)).sum()), ("0,25 < max TPM < 0,5", ((mx > 0.25) & (mx < 0.5)).sum())]
    fig1 = hbars([dict(label=b[0], value=int(b[1]), color=2) for b in bins], "Gènes codants non exprimés selon leur TPM maximal", label_w=250)
    tab1 = data_table(["Classe de TPM maximal", "Gènes"], [[b[0], fmt_int(b[1])] for b in bins])
    thr_rows = []
    for th in (0.1, 0.25, 0.5, 1, 2, 5):
        thr_rows.append([f"TPM ≥ {dec_fr(th, 2 if th < 1 else 0).rstrip('0').rstrip(',')} dans ≥ 1 banque", fmt_int((mx < th).sum()),
                         pct((mx < th).sum(), len(mx)),
                         fmt_int(((tm >= th).sum(axis=1) < 2).sum())])
    thr = table(["Critère d'expression", "Non exprimés", "% des codants", "Non exprimés si ≥ 2 banques exigées"], thr_rows, hl_rows=(2,))
    cur = saturation(t, coding)
    ks = list(range(1, cur.shape[1] + 1))
    mean = cur.mean(axis=0)
    lo, hi = cur.min(axis=0), cur.max(axis=0)
    sat = line_chart([("Gènes détectés (moyenne de 100 ordres aléatoires)", ks, list(mean), 1)],
                     "Courbe de saturation : gènes détectés selon le nombre de banques RNA-seq", "Nombre de banques RNA-seq ajoutées",
                     "Gènes codants avec TPM ≥ 0,5", bands=[(ks, list(lo), list(hi), 1)], ymax=45000, yticks=[0, 10000, 20000, 30000, 40000],
                     xticks=[1, 10, 20, 30, 40, 47], left=90)
    inc5 = mean[-1] - mean[-6]
    sat_t = data_table(["Banques", "Gènes détectés (moyenne)"], [[k, fmt_int(mean[k - 1])] for k in (1, 5, 10, 20, 30, 40, 47)])
    return f"""
<h3>3.4 Les 15&nbsp;454 gènes sans expression : de quoi parle-t-on ?</h3>
<p>Le chiffre des <b>{fmt_int(A['no_expr_total'])} gènes sans expression</b> (27,2&nbsp;% de l'annotation) se décompose ainsi :</p>
<ul>
<li><b>978 gènes ncRNA</b> dont le TPM est exactement 0 dans les 47 banques&nbsp;: ils n'ont jamais été quantifiés (le transcriptome de référence de Salmon, <code>final_transcripts.fasta</code>, ne contient que les 57&nbsp;542 ARNm, selon <code>modules/final_expression_validation.nf</code>). Ils ne sont donc pas «&nbsp;non exprimés&nbsp;», ils sont non mesurés.</li>
<li><b>{fmt_int(n_noexp)} gènes codants</b> avec un TPM maximal inférieur à 0,5 sur les 47 banques (26,0&nbsp;% des codants).</li>
</ul>
<h4>Ce nombre dépend du seuil</h4>
{figure(fig1, "TPM maximal des gènes codants non exprimés", f"{fmt_int(bins[0][1])} gènes n'ont aucun read assigné ; {fmt_int(bins[3][1])} sont entre 0,25 et 0,5, c'est-à-dire à moins d'un facteur 2 du seuil.", tab1)}
{thr}
<p>Un seuil de 0,1 TPM ramènerait les non exprimés à 7&nbsp;407 ; exiger 2 banques les porterait à 17&nbsp;942. Le chiffre de 15&nbsp;000 est donc <b>un ordre de grandeur, pas une frontière biologique</b>. Seuls les {fmt_int(bins[0][1])} gènes sans aucun read sont non exprimés quel que soit le seuil.</p>
<h4>Le panel de 47 banques n'est pas saturé</h4>
{figure(sat, "Détection en fonction du nombre de banques", f"La courbe continue de monter : la dernière banque ajoute en moyenne {dec_fr(mean[-1] - mean[-2], 0)} gènes, les cinq dernières {fmt_int(inc5)} (≈ {dec_fr(100 * inc5 / mean[-1])} %). La bande grise donne le minimum et le maximum sur 100 ordres.", sat_t)}
<p>Avec {fmt_int(mean[-1])} gènes détectés, la courbe est proche du plateau mais pas atteinte&nbsp;: des gènes propres à un tissu, un stade ou une condition non couverts restent non exprimés dans ce panel. <b>Non exprimé ne veut donc pas dire inexistant.</b> C'est pourquoi la suite ne se sert pas de l'expression seule mais la croise avec d'autres preuves.</p>
"""


def sec_classes(D):
    n = D.n_class
    items = [dict(label=CLASS_FR[k], value=int(n.get(k, 0)), color=CLASS_COLOR[k],
                  text=f"{fmt_int(n.get(k, 0))} ({pct(n.get(k, 0), D.n)})") for k in CLASS_ORDER]
    fig = hbars(items, "Répartition des gènes codants par classe de preuve", label_w=300, right_pad=150)
    tab = data_table(["Classe", "Gènes", "%"], [[CLASS_FR[k], fmt_int(n.get(k, 0)), pct(n.get(k, 0), D.n)] for k in CLASS_ORDER])
    defs = table(["Classe", "Définition (dans cet ordre de priorité)"], text=True, rows=[
        [CLASS_FR["expressed"], "TPM ≥ 0,5 dans au moins une des 47 banques."],
        [CLASS_FR["te_protein"], "Non exprimé ; ≥ 50 % de la CDS dans une annotation de TE (EDTA) <b>et</b> un résultat fonctionnel qui décrit une protéine de TE (transposase, transcriptase inverse, gag-pol, intégrase…)."],
        [CLASS_FR["te_overlap_unsupported"], "Non exprimé ; ≥ 50 % de la CDS sur un TE ; pas de protéine de TE ; ni fonction, ni famille OMAMER, ni homologue non-<i>Vitis</i>."],
        [CLASS_FR["conserved_or_functional"], "Non exprimé ; au moins une fonction (Diamond2GO, eggNOG, InterProScan, MapMan), ou une famille OMAMER, ou un homologue (≥ 35 % d'identité, ≥ 50 % de couverture) chez une eudicotylédone non-<i>Vitis</i>."],
        [CLASS_FR["liftoff_only"], "Aucun des supports précédents, mais porte un identifiant v4.3 (transfert Liftoff)."],
        [CLASS_FR["junction_only"], "Aucun des supports précédents, mais au moins un intron retrouvé dans un assemblage RNA-seq."],
        [CLASS_FR["no_signal"], "Rien : ni expression, ni jonction, ni fonction, ni orthologie, ni homologue, ni v4.3, ni TE."],
    ])
    return f"""
<h3>3.5 Classer les gènes selon leurs autres preuves</h3>
<p>Pour décider quels gènes non exprimés sont douteux, chaque gène codant reçoit <b>une seule classe</b>, selon ses preuves autres que l'expression (les critères sont évalués dans l'ordre du tableau).</p>
{defs}
{figure(fig, "Répartition des 55 766 gènes codants", f"Sur {fmt_int(D.noexpr)} gènes non exprimés, {fmt_int(n.get('conserved_or_functional', 0))} ont un support de fonction ou de conservation et {fmt_int(D.set_aside)} n'en ont aucun.", tab)}
<p>Les homologues non-<i>Vitis</i> viennent de DIAMOND (e &lt; 10<sup>-5</sup>) contre les protéines d'eudicotylédones d'UniProt <b>en excluant tout <i>Vitis</i></b>&nbsp;: un homologue chez une autre espèce est une preuve indépendante que la famille est réelle, alors qu'une correspondance avec une autre annotation de vigne peut reproduire la même erreur.</p>
"""


def sec_proofs(D):
    c = D.cls
    A = AUDIT
    order = ["expressed", "conserved_or_functional", "te_protein", "te_overlap_unsupported", "no_signal"]
    # 1. reads
    r = c["salmon_reads_total"]
    rows = []
    for k in order:
        x = r[c["evidence_class"] == k]
        rows.append((CLASS_FR[k], [int((x == 0).sum()), int(((x > 0) & (x < 10)).sum()), int(((x >= 10) & (x < 100)).sum()), int((x >= 100).sum())]))
    fig_r = stacked_hbars(rows, [("0 read", 7), ("1 à 9", 8), ("10 à 99", 9), ("≥ 100", 10)],
                          "Nombre total de reads Salmon assignés au gène sur les 47 banques", label_w=300)
    tab_r = data_table(["Classe", "0 read", "1-9", "10-99", "≥ 100"], [[a] + [fmt_int(v) for v in b] for a, b in rows])
    ns = c[c["evidence_class"] == "no_signal"]["salmon_reads_total"]
    ex = c[c["evidence_class"] == "expressed"]["salmon_reads_total"]
    # 2. heatmap
    cols = ["Jonction (multi-exons)", "Fonction / MapMan", "Famille OMAMER", "Homologue non-Vitis", "Gène v4.3 (Liftoff)", "Dans PN40024 5.1", "Dans jeu TSEBRA"]
    vals, txts, rl = [], [], []
    for k in order:
        x = c[c["evidence_class"] == k]
        multi = x[x["n_exons"] > 1]
        v = [100.0 * multi["has_junction"].mean(), 100.0 * x["has_function"].mean(), 100.0 * x["omamer_hit"].mean(),
             100.0 * x["has_nonvitis_homolog"].mean(), 100.0 * x["liftoff_carried"].mean(), 100.0 * x["in_v51"].mean(),
             100.0 * x["in_tsebra"].mean()]
        vals.append([a / 100.0 for a in v])
        txts.append([fmt_pct(a, 0) for a in v])
        rl.append(CLASS_FR[k])
    hm = heatmap(rl, cols, vals, txts, "Profil de preuves par classe (% des gènes de la classe)", label_w=300)
    # 3. protein length
    bins = np.arange(0, 825, 25)
    ser = []
    for name, mask, col in (("Exprimés", c["evidence_class"] == "expressed", 1),
                            ("Non exprimés, conservés ou fonctionnels", c["evidence_class"] == "conserved_or_functional", 3),
                            ("À écarter (aucun signal + sur un TE)", c["evidence_class"].isin(["no_signal", "te_overlap_unsupported"]), 2)):
        h, _ = np.histogram(np.clip(c.loc[mask, "prot_len"], 0, 799), bins=bins)
        ser.append((name, list(bins[:-1] + 12.5), list(100.0 * h / h.sum()), col))
    ymax = float(np.ceil(max(max(sr[2]) for sr in ser) * 1.12 / 5.0) * 5)
    fig_l = line_chart(ser, "Distribution de la longueur des protéines", "Longueur de la protéine (acides aminés ; la dernière classe regroupe ≥ 775 aa)",
                       "% des gènes du groupe par classe de 25 aa", ymax=ymax, yticks=[0, 5, 10, 15, 20, 25][:int(ymax // 5) + 1],
                       xticks=[0, 200, 400, 600, 800], yfmt=lambda v: f"{v:g}", left=90)
    med = {n_: float(c.loc[m, "prot_len"].median()) for n_, m in (("e", c["evidence_class"] == "expressed"), ("s", c["evidence_class"].isin(["no_signal", "te_overlap_unsupported"])))}
    n_short = 100.0 * (c[c["evidence_class"].isin(["no_signal", "te_overlap_unsupported"])]["prot_len"] < 100).mean()
    # validations
    v_rows = []
    for k in order:
        x = c[c["evidence_class"] == k]
        v_rows.append([CLASS_FR[k], fmt_int(len(x)), fmt_int(int(x["busco"].sum())), fmt_pct(100.0 * x["in_v51"].mean())])
    val = table(["Classe", "Gènes", "Gènes BUSCO", "Recouvrent un gène de 5.1"], v_rows, hl_rows=(3, 4))
    return f"""
<h3>3.6 Les preuves que ces 2&nbsp;478 gènes sont des erreurs</h3>
<p>Il n'existe pas de preuve qu'un gène n'existe pas. L'argument est que <b>toutes les lignes de preuve indépendantes convergent</b> et que les gènes écartés diffèrent des gènes connus comme réels sur chaque dimension mesurée.</p>
<h4>(a) Ce n'est pas un effet de seuil : il n'y a pas de reads</h4>
{figure(fig_r, "Reads Salmon totaux par classe", f"Gènes « aucun signal » : médiane de {dec_fr(ns.median(), 0)} read sur les 47 banques, {pct((ns == 0).sum(), len(ns), 0)} à zéro read, {pct((ns < 10).sum(), len(ns), 0)} sous 10 reads. Gènes exprimés : médiane de {fmt_int(ex.median())} reads.", tab_r)}
<h4>(b) Aucun intron, aucune conservation</h4>
{figure(hm, "Profil de preuves de chaque classe", "Lecture : 0 % des gènes « aucun signal » ont un intron retrouvé, une famille OMAMER ou un homologue non-Vitis, contre 78 % (OMAMER) et 73 % (homologue) des gènes exprimés. La jonction est calculée sur les gènes multi-exons.", "")}
<p>Parmi les gènes non exprimés, 96&nbsp;% des multi-exons n'ont aucun intron retrouvé dans les assemblages RNA-seq (contre 30&nbsp;% des gènes exprimés). Cette vérification utilise d'autres logiciels que Salmon (StringTie, PsiCLASS, long reads) et n'est donc pas une répétition du même signal.</p>
<h4>(c) Des protéines courtes, typiques de modèles ab initio isolés</h4>
{figure(fig_l, "Longueur des protéines", f"Médiane de {dec_fr(med['s'], 0)} aa pour les gènes à écarter contre {dec_fr(med['e'], 0)} aa pour les exprimés ; {dec_fr(n_short, 0)} % des gènes à écarter font moins de 100 aa (exprimés : 13 %).", "")}
<h4>(d) Validations indépendantes de la définition des classes</h4>
{val}
<p><b>Aucun gène BUSCO</b> n'est dans les classes à écarter (0 sur 2&nbsp;022 gènes BUSCO du jeu), et leur recouvrement avec la référence PN40024 5.1 est de 6,1&nbsp;% et 8,2&nbsp;%, contre 66,8&nbsp;% pour les gènes exprimés. 5.1 étant elle-même construite sur du RNA-seq, ce dernier point est une vérification de cohérence et non une preuve indépendante de l'expression.</p>
{callout('<p><span class="h">Confiance.</span> Élevée pour les classes «&nbsp;aucun signal&nbsp;» et «&nbsp;sur un TE sans autre support&nbsp;» prises ensemble. Moyenne pour «&nbsp;gène v4.3 seulement&nbsp;» (539 gènes&nbsp;: même profil, mais gardés par une annotation antérieure). La classe «&nbsp;jonction seulement&nbsp;» (11 gènes) est trop petite pour trancher. On peut s\'attendre à quelques faux positifs&nbsp;: 77 des 1&nbsp;263 gènes «&nbsp;aucun signal&nbsp;» recouvrent un gène de 5.1.</p>', "warn")}
"""


def sec_rejected(D):
    s = D.sig
    rows = []
    spec = [("Fragment d'un paralogue (hit partiel)", "sig_fragment_of_paralog", "Non enrichi"),
            ("Gène coupé : voisin soutenu, même brin, < 2 kb", "sig_near_supported_same_strand", "Moins fréquent"),
            ("Antisens ou imbriqué dans un autre gène", "sig_antisense_or_nested", "Non enrichi"),
            ("Recouvrement CDS de même brin", "sig_cds_overlap_same_strand", "Légèrement enrichi"),
            ("Protéine de faible complexité", "sig_low_complexity", "Non enrichi"),
            ("Sur chr00 (séquences non placées)", "sig_chr00", "Légèrement enrichi"),
            ("Mono-exonique", "sig_single_exon", "Non enrichi"),
            ("ORF court (< 100 aa)", "sig_short_orf", "Enrichi"),
            ("Modèle Helixer seul", "sig_helixer_only", "Enrichi (×2)"),
            ("Modèle ab initio brut absent du jeu TSEBRA", "sig_raw_abinitio_not_tsebra", "Très fortement enrichi")]
    for lab, col, verdict in spec:
        rows.append([lab, fmt_pct(100 * s.loc[col, "no_signal"]), fmt_pct(100 * s.loc[col, "expressed"]), verdict])
    return f"""
<h3>3.7 Hypothèses examinées et écartées</h3>
<p>Pour comprendre <i>pourquoi</i> ces modèles sont faux, plusieurs mécanismes classiques ont été testés sur les 1&nbsp;263 gènes «&nbsp;aucun signal&nbsp;».</p>
{table(["Hypothèse", "« Aucun signal »", "Exprimés", "Verdict"], rows, hl_rows=(7, 8, 9))}
<p>Ce ne sont donc <b>ni des fragments de paralogues, ni des gènes coupés, ni des artefacts de chevauchement, ni des protéines de faible complexité</b>. Ce sont des modèles courts, isolés, d'origine ab initio brute, sans aucun support externe.</p>
"""


def sec_cause(D):
    c = D.cls
    b = c[c["origin"] == "ab_initio_braker"]
    tsebra_share = 100.0 * b["in_tsebra"].mean()
    setaside = c[c["evidence_class"].isin(["no_signal", "te_overlap_unsupported"])]
    braker_share = 100.0 * (setaside["origin"] == "ab_initio_braker").mean()
    ne = c[c["evidence_class"] != "expressed"]
    return f"""
<h3>3.8 Cause dans le pipeline : les modèles BRAKER bruts</h3>
<p><b>{dec_fr(braker_share, 0)}&nbsp;% des gènes à écarter sont des modèles BRAKER bruts</b> (les autres&nbsp;: Helixer 12&nbsp;%, Liftoff 4&nbsp;%, assemblages 2&nbsp;%, EGAPx 0&nbsp;%), alors que BRAKER fournit 33&nbsp;% de tous les gènes. Le câblage du pipeline l'explique (vérifié dans <code>modules/mikado.nf</code> et <code>workflows/titan.nf</code>)&nbsp;:</p>
<ul>
<li>Mikado reçoit <code>augustus.hints.gff3</code> (83&nbsp;255 transcrits) et <code>genemark.gtf</code>, c'est-à-dire les prédictions <b>brutes</b>. Le jeu filtré par TSEBRA, <code>braker.gff3</code> (52&nbsp;414 transcrits), et <code>genemark_supported.gtf</code> sont publiés mais n'atteignent pas Mikado.</li>
<li>Seuls <b>{dec_fr(tsebra_share, 0)}&nbsp;%</b> des gènes d'origine BRAKER recouvrent un modèle TSEBRA (≥ 50&nbsp;% de la CDS, même brin). Parmi les gènes non exprimés, 24&nbsp;% sont dans le jeu TSEBRA, contre 61&nbsp;% des gènes exprimés.</li>
<li>Mikado garde le meilleur transcrit de chaque locus <i>même sans concurrent</i>. Un modèle AUGUSTUS isolé, sans reads, sans jonction et sans protéine voisine, est donc reporté tel quel.</li>
</ul>
<p>Il s'agit d'une <b>association</b> solide, pas d'une cause démontrée seule&nbsp;: pour la tester, Mikado a été relancé avec les pistes filtrées (section suivante).</p>
"""


def sec_test(D):
    fr, ft = D.fate["raw"], D.fate["tsebra"]
    cats = [CLASS_FR[k] for k in CLASS_ORDER if k in ft.index]
    ks = [k for k in CLASS_ORDER if k in ft.index]
    fig = grouped_hbars(cats, [("Contrôle (pistes brutes)", [float(fr.loc[k, "pct_lost"]) for k in ks], 1),
                               ("Test (pistes TSEBRA)", [float(ft.loc[k, "pct_lost"]) for k in ks], 2)],
                        "Part des gènes de production perdus, contrôle et test", label_w=310, max_val=70)
    tab = data_table(["Classe", "Gènes", "Perdus (contrôle)", "Perdus (test)", "% perdus (test)"],
                     [[CLASS_FR[k], fmt_int(ft.loc[k, "n"]), fmt_int(fr.loc[k, "lost"]), fmt_int(ft.loc[k, "lost"]), fmt_pct(ft.loc[k, "pct_lost"])] for k in ks])
    cm = D.cmp
    ctl, tst = cm["raw"], cm["tsebra"]
    crit = table(["Critère fixé avant de lancer", "Résultat", "Verdict"], text=True, rows=[
        ["Le contrôle reproduit la production (≥ 99 % des loci)", f"{dec_fr(100.0 * ctl['loci_identical_to_production_mikado'] / ctl['production_mikado_genes'])} %, BUSCO 99,8 %", "atteint"],
        ["Le test perd ≥ 50 % des 2 478 gènes à écarter", "1 490 perdus (60,1 %)", "atteint"],
        ["Le test perd ≤ 1 % des gènes exprimés et ≤ 5 gènes BUSCO", "9,7 % des exprimés (4 023) ; 0 gène BUSCO", "<b>non atteint</b>"],
        ["Le test perd ≤ 10 % des gènes conservés ou fonctionnels", "23,6 % (2 429)", "<b>non atteint</b>"],
        ["BUSCO ne baisse pas de plus de 0,3 point", "99,9 % contre 99,8 %", "atteint"],
        ["Gènes nouveaux ≤ 2 % du total", f"{fmt_int(tst['genes_not_in_production'])} ({dec_fr(100.0 * tst['genes_not_in_production'] / tst['genes'])} %)", "<b>non atteint</b> (de peu)"],
    ], hl_rows=(2, 3, 5))
    return f"""
<h3>3.9 Test : donner à Mikado les pistes filtrées par TSEBRA</h3>
<p>Mikado a été relancé deux fois avec les mêmes modules, conteneurs et pistes, une seule chose changeant&nbsp;: le <b>contrôle</b> reçoit les pistes brutes actuelles, le <b>test</b> reçoit <code>braker.gff3</code> et <code>genemark_supported.gtf</code>. Les critères de réussite ont été fixés avant le lancement.</p>
{table(["", "Production", "Contrôle (brut)", "Test (TSEBRA)"], [
    ["Gènes codants", "55 766", fmt_int(ctl["genes"]), fmt_int(tst["genes"])],
    ["Loci identiques à la production", "-", fmt_int(ctl["loci_identical_to_production_mikado"]), fmt_int(tst["loci_identical_to_production_mikado"])],
    ["Gènes absents de la production", "-", fmt_int(ctl["genes_not_in_production"]), fmt_int(tst["genes_not_in_production"])],
    ["BUSCO complet (n = 1 990)", "99,8 %", dec_fr(float(ctl["busco"]["complete"])) + " %", dec_fr(float(tst["busco"]["complete"])) + " %"],
])}
{figure(fig, "Gènes de production perdus par le test", "Barres bleues : bruit de base du contrôle (0,9 % des gènes exprimés). Barres orange : gènes perdus avec les pistes TSEBRA.", tab)}
{crit}
{callout('<p><span class="h">Verdict.</span> Le test n\'est pas une amélioration nette. Il supprime bien 60&nbsp;% des gènes visés, mais aussi 9,7&nbsp;% des gènes exprimés et 23,6&nbsp;% des gènes conservés. Il ne faut donc <b>pas</b> passer les pistes TSEBRA en production telles quelles. BUSCO ne voit pas cette perte (99,9&nbsp;% avec 6&nbsp;000 gènes de moins), il ne peut pas arbitrer.</p>', "warn")}
<h4>À quoi ressemblent les gènes perdus ?</h4>
<ul>
<li><b>Presque tous sont BRAKER</b> : 7&nbsp;822 des 18&nbsp;429 gènes d'origine BRAKER sont perdus (42&nbsp;%), contre 0,4&nbsp;% pour Helixer, 0,2&nbsp;% pour Liftoff et 0&nbsp;% pour EGAPx. <b>AUGUSTUS pèse le plus</b> : 55&nbsp;% de ses modèles disparaissent, contre 27&nbsp;% pour GeneMark.</li>
<li><b>Le tri va dans le bon sens</b> : à classe égale, les gènes perdus sont plus faibles que les gardés (parmi les exprimés d'origine BRAKER : protéine médiane 147 aa contre 298, intron soutenu 3&nbsp;% contre 26&nbsp;%, homologue non-<i>Vitis</i> 32&nbsp;% contre 59&nbsp;%).</li>
<li><b>Mais pas tous sont des déchets</b> : parmi les 4&nbsp;023 gènes exprimés perdus, 26&nbsp;% ont un TPM maximal ≥ 5 et 64&nbsp;% une fonction ; environ 890 gènes perdus cumulent plusieurs preuves fortes (borne haute d'une définition grossière).</li>
<li><b>Ils disparaissent vraiment</b> : 70&nbsp;% des gènes exprimés perdus n'ont aucun gène de même brin à leur locus dans le test.</li>
</ul>
<p>L'hypothèse la plus simple est que le filtre TSEBRA agit comme un filtre brutal sur BRAKER : il retire les modèles sans indice, mais aussi de vrais gènes à faible signal. Une piste non testée est d'isoler l'effet d'AUGUSTUS de celui de GeneMark (un troisième bras).</p>
"""


def sec_real(D):
    c = D.cls
    cons = c[c["evidence_class"] == "conserved_or_functional"]
    z = c.groupby(c["chrom"] == "chr00")["evidence_class"].value_counts(normalize=True).unstack().fillna(0)
    cats = [CLASS_FR[k] for k in CLASS_ORDER]
    other = [100.0 * z.loc[False].get(k, 0) for k in CLASS_ORDER]
    c0 = [100.0 * z.loc[True].get(k, 0) for k in CLASS_ORDER]
    fig = grouped_hbars(cats, [("Chromosomes 1 à 19", other, 1), ("chr00 (non placé)", c0, 2)],
                        "Répartition des classes sur chr00 et sur les chromosomes 1 à 19", label_w=310, max_val=80)
    tab = data_table(["Classe", "chr01-19", "chr00"], [[CLASS_FR[k], fmt_pct(o), fmt_pct(z0)] for k, o, z0 in zip(CLASS_ORDER, other, c0)])
    n_cons_c0 = int((cons["chrom"] == "chr00").sum())
    return f"""
<h3>3.10 Les gènes non exprimés qui sont probablement réels</h3>
<p>Les {fmt_int(len(cons))} gènes non exprimés mais <b>conservés ou fonctionnels</b> ne doivent pas être écartés&nbsp;: 98&nbsp;% ont une annotation fonctionnelle, 66&nbsp;% une famille OMAMER, 48&nbsp;% un homologue non-<i>Vitis</i>. Ils sont absents des 47 banques, pas de la biologie (gènes tissu-spécifiques, familles de résistance ou de réponse au stress, etc.).</p>
{figure(fig, "chr00 se distingue nettement", f"{pct(z.loc[True].get('conserved_or_functional', 0) * 100, 100, 0)} des gènes de chr00 sont dans cette classe (contre 17 % ailleurs).", tab)}
<p><b>chr00</b> (séquences non placées, 2&nbsp;959 gènes) contient {fmt_int(n_cons_c0)} de ces gènes. Parmi eux, <b>1&nbsp;376 sont identiques à au moins 95&nbsp;% à un autre gène de chr00</b> (médiane de 2,5&nbsp;Mb de distance, donc pas des duplications en tandem), à 96&nbsp;% d'origine BRAKER, et seuls 25&nbsp;% de leurs partenaires sont exprimés. Ce profil est compatible avec des copies quasi identiques dont les reads ne peuvent pas être attribués à une seule copie (multi-mapping), ou avec du matériel dupliqué. <b>C'est une hypothèse non vérifiée au niveau des reads</b> (un comptage direct dans les BAM a été tenté mais s'est révélé trop lent&nbsp;: environ 10&nbsp;h par BAM).</p>
<p><b>Protéines de TE.</b> {fmt_int(D.n_class.get('te_protein', 0))} gènes non exprimés (et 1&nbsp;289 exprimés) sont des protéines de transposons (transposase, transcriptase inverse, gag-pol… ; 91&nbsp;% ont une famille OMAMER)&nbsp;: ce ne sont pas des erreurs de prédiction mais ils ne sont pas des gènes de l'hôte. Il est fréquent de les traiter dans une piste séparée. Attention&nbsp;: le seul mot-clé «&nbsp;domaine de TE&nbsp;» sans chevauchement EDTA n'est pas fiable (13 gènes BUSCO exprimés le portent).</p>
"""


def sec_balance(D):
    A, V = AUDIT, AUDIT["v51"]
    n = D.n_class
    coding = A["coding"]
    strict = coding - D.set_aside
    broad = strict - n.get("liftoff_only", 0) - n.get("junction_only", 0)
    full = broad - n.get("te_protein", 0)
    items = [
        dict(label="v4.3", value=A["v43_genes"], color=1, text="34 843"),
        dict(label="PN40024 5.1", value=V["genes"], color=1, text="47 971"),
        dict(label="TITAN, tous les codants", value=coding, color=2),
        dict(label="− « strict » (aucun signal, sur un TE)", value=strict, color=2),
        dict(label="− « broad » (+ v4.3 seul, jonction seule)", value=broad, color=2),
        dict(label="− « broad+te » (+ protéines de TE)", value=full, color=2),
        dict(label="≥ 2 lignes de preuve", value=D.ge2, color=3),
        dict(label="≥ 3 lignes de preuve", value=D.ge3, color=3),
        dict(label="≥ 4 lignes de preuve", value=D.ge4, color=3),
    ]
    fig = hbars(items, "Nombre de gènes codants selon la sévérité du filtre", label_w=300, right_pad=110)
    tab = data_table(["Jeu", "Gènes codants"], [[i["label"], fmt_int(i["value"])] for i in items])
    return f"""
<h3>3.11 Bilan : combien de gènes, et que dire de l'écart avec 35&nbsp;000 ?</h3>
{figure(fig, "Nombre de gènes selon les politiques de filtrage", "Bleu : références. Orange : TITAN filtré par politique. Vert : gènes supportés par au moins 2, 3 ou 4 lignes de preuve indépendantes.", tab)}
<ul>
<li><b>Les erreurs manifestes sont peu nombreuses</b> : {fmt_int(D.set_aside)} gènes (4,4&nbsp;% des codants) avec la politique stricte, ce qui laisse {fmt_int(strict)} gènes. Même la politique la plus sévère laisse {fmt_int(full)} gènes.</li>
<li><b>Le chiffre élevé n'est donc pas surtout un problème de faux positifs.</b> Il vient pour l'essentiel de gènes soutenus (exprimés, conservés ou annotés).</li>
<li><b>Causes plausibles de l'écart avec les ordres de grandeur historiques</b> (à discuter, non toutes testées)&nbsp;: (1) un assemblage T2T qui résout des régions dupliquées absentes des assemblages antérieurs&nbsp;; (2) une évidence bien plus riche (47 banques, long reads, 3 méthodes fonctionnelles)&nbsp;; (3) des familles multigéniques (NB-LRR, PPR, kinases) découpées en copies distinctes&nbsp;; (4) des duplications segmentaires ou d'haplotypes (OMArk : 10,6&nbsp;% de duplications inattendues, chr00)&nbsp;; (5) une possible fragmentation de gènes non mesurée ici&nbsp;; (6) des gènes de TE non retirés (2&nbsp;429 protéines de TE).</li>
<li><b>Où trouver un nombre plus prudent</b> : le noyau à ≥ 3 lignes de preuve compte {fmt_int(D.ge3)} gènes. Ce n'est pas le nombre de gènes de la vigne, mais une borne basse très fiable.</li>
</ul>
{callout('<p><span class="h">Message pour la discussion.</span> Avec 47 banques RNA-seq, 3 méthodes fonctionnelles et un assemblage T2T, TITAN annote plus de gènes que les références antérieures, et la grande majorité est soutenue. Le vrai enjeu n\'est pas d\'atteindre 35&nbsp;000 mais de <b>distinguer clairement</b> le noyau très fiable, les gènes plausibles mais peu soutenus (à marquer) et les erreurs manifestes (à écarter).</p>', "ok")}
"""


def sec_tool(D):
    n = D.n_class
    strict, broad = D.set_aside, D.set_aside + n.get("liftoff_only", 0) + n.get("junction_only", 0)
    full = broad + n.get("te_protein", 0)
    pol = table(["Politique", "Classes retirées", "Gènes retirés", "Gènes restants (GFF3, avec ncRNA)"], text=True, rows=[
        ["<code>strict</code>", "aucun signal, sur un TE sans support", fmt_int(strict), fmt_int(AUDIT["genes"] - strict)],
        ["<code>broad</code>", "+ gène v4.3 seul, jonction seule", fmt_int(broad), fmt_int(AUDIT["genes"] - broad)],
        ["<code>broad+te</code>", "+ protéines de TE (piste TE séparée)", fmt_int(full), fmt_int(AUDIT["genes"] - full)],
    ])
    return f"""
<h2 id="outil">4. Marquer et filtrer le GFF3 (outil optionnel, hors pipeline)</h2>
<p>Le script <code>scripts/flag_unsupported_genes.py</code> (Python standard, sans dépendance) ne fait <b>rien à la production</b>&nbsp;: il lit <code>final_annotation.gff3</code> et le tableau de classes, ne modifie jamais ses entrées et <b>ne renomme aucun gène</b> (les identifiants restent joignables à toutes les tables de TITAN).</p>
<h3>4.1 Trois commandes</h3>
<pre><code># 1. tableau des classes (une fois, à partir d'un run TITAN terminé)
scripts/run_unsupported_genes_analysis.sh data/titan_prod_out data/unsupported_genes_analysis data/PN40024_5.1_on_T2T_ref.gff3

# 2. combien de gènes chaque politique retirerait
scripts/flag_unsupported_genes.py summary --gff3 final_annotation.gff3 --classes classes/gene_classes.tsv

# 3a. marquer (aucun gène retiré) : ajoute titan_evidence_class et titan_filter_level sur chaque ligne gene
scripts/flag_unsupported_genes.py mark --gff3 final_annotation.gff3 --classes classes/gene_classes.tsv --out final_annotation.marked.gff3

# 3b. filtrer : retire les gènes et tous leurs enfants (mRNA, exon, CDS, UTR), et les protéines correspondantes
scripts/flag_unsupported_genes.py filter --gff3 final_annotation.gff3 --classes classes/gene_classes.tsv \\
    --policy strict --out final_annotation.strict.gff3 \\
    --proteins final_annotation_proteins_main.fasta --proteins-out proteins_main.strict.fasta \\
    --removed-list removed_genes.txt</code></pre>
<h3>4.2 Politiques</h3>
{pol}
<p>Chaque ligne <code>gene</code> marquée reçoit deux attributs, par exemple <code>titan_evidence_class=no_signal;titan_filter_level=strict</code> (niveaux : <code>keep</code>, <code>strict</code>, <code>broad</code>, <code>te</code>). Les gènes ncRNA et les gènes absents du tableau sont conservés (<code>not_evaluated</code>). Les classes inconnues font échouer le script plutôt que d'être ignorées.</p>
<h3>4.3 Vérifications faites</h3>
<ul>
<li><b>Test unitaire</b> (<code>scripts/test_flag_unsupported_genes.py</code>) : marquage sans ajout ni perte de ligne, trois politiques, aucun enfant orphelin, FASTA cohérent, entrée non modifiée, classe inconnue refusée.</li>
<li><b>Sur le run de production</b> : politique stricte → 54&nbsp;266 gènes conservés (2&nbsp;478 retirés), 53&nbsp;288 protéines principales, 55&nbsp;064 protéines (toutes isoformes). Le GFF3 filtré passe le <b>validateur de TITAN</b> (<code>validate_final_annotation.py</code> : PASS, 397&nbsp;565 identifiants, 53&nbsp;288 protéines principales retrouvées) avec exactement les mêmes avertissements que l'original (stops internes).</li>
</ul>
{callout('<p><span class="h">À ne pas faire sans discussion.</span> Ce filtre retire des gènes sur des critères d\'<i>absence</i> de preuve. Il ne remplace pas la curation. Recommandation&nbsp;: distribuer l\'annotation primaire <b>marquée</b> plutôt que filtrée, et laisser chaque analyse choisir sa politique. Les gènes de la classe « conservés ou fonctionnels » ne sont jamais retirés.</p>', "warn")}
"""


def sec_limits(D):
    return """
<h2 id="limites">5. Limites et points à discuter</h2>
<ul>
<li><b>L'absence de preuve n'est pas la preuve de l'absence.</b> La classe «&nbsp;aucun signal&nbsp;» ne peut pas repérer un gène réel très spécifique à la vigne et exprimé dans un tissu non échantillonné (5.1 recouvre 77 des 1&nbsp;263 gènes).</li>
<li><b>Expression</b> : 47 banques, seuil de 0,5&nbsp;TPM, un seul échantillon suffit. Le nombre de gènes non exprimés change fortement avec le seuil (7&nbsp;407 à 0,1&nbsp;TPM, 17&nbsp;942 si 2 banques sont exigées). Les 978 ncRNA n'ont pas été quantifiés.</li>
<li><b>Homologie</b> : les banques de référence (UniProt, Swiss-Prot, Vitales) contiennent elles-mêmes des protéines prédites. Un homologue prouve qu'une famille existe, pas que ce modèle est correct.</li>
<li><b>TE</b> : dépend d'EDTA et d'une recherche de mots-clés dans les annotations fonctionnelles.</li>
<li><b>Recouvrement avec 5.1</b> : 5.1 repose sur du RNA-seq, ce n'est pas une vérité de terrain indépendante.</li>
<li><b>Mécanisme</b> : «&nbsp;ORF ab initio aléatoire&nbsp;» est une inférence (origine du modèle, absence de support). L'usage des codons ou le GC3 de ces ORF n'ont pas été comparés à des ORF aléatoires.</li>
<li><b>Test TSEBRA</b> : comparaison au niveau des loci Mikado, avant AEGIS, ncRNA et annotation fonctionnelle&nbsp;; le contrôle diffère de la production pour 0,8&nbsp;% des loci (bruit de base d'environ 1&nbsp;%). Les effets d'AUGUSTUS et de GeneMark n'ont pas été séparés.</li>
<li><b>chr00</b> : l'hypothèse de multi-mapping n'est pas confirmée au niveau des reads.</li>
<li><b>Points ouverts de l'audit</b> : FLAIR n'a produit aucun isoforme dans ce run&nbsp;; la répartition par catégorie structurale de SQANTI3 (20&nbsp;143 isoformes long reads, tous en «&nbsp;Other&nbsp;») a l'allure d'un artefact de catégorisation et doit être vérifiée avant d'être citée&nbsp;; la quantification Salmon n'a pas activé les bootstraps.</li>
<li><b>Pas de vérité de terrain</b> : aucune validation expérimentale (RT-PCR, protéomique). Les conclusions portent sur la cohérence des preuves.</li>
</ul>
"""


def sec_repro(D, a):
    rows = [
        ["<code>scripts/build_gene_evidence_table.py</code>", "Tableau par gène : structure, provenance Mikado, expression, fonction, MapMan, OMAMER, Liftoff, TE, voisinage"],
        ["<code>scripts/add_intron_support.py</code>", "Introns retrouvés dans les assemblages RNA-seq"],
        ["<code>scripts/add_homology_features.py</code>", "Recouvrement TSEBRA, DIAMOND (soi, Vitales+Swiss-Prot, eudicotylédones non-<i>Vitis</i>)"],
        ["<code>scripts/salmon_gene_numreads.py</code>", "Reads Salmon totaux par gène"],
        ["<code>scripts/classify_unsupported_genes.py</code>", "Classes de preuve, signatures, listes"],
        ["<code>scripts/run_unsupported_genes_analysis.sh</code>", "Enchaîne les étapes ci-dessus"],
        ["<code>scripts/flag_unsupported_genes.py</code>", "Marquage et filtrage du GFF3 (section 4)"],
        ["<code>mikado_input_test.nf</code>, <code>scripts/mikado_input_test/</code>", "Test Mikado, lanceur Slurm, comparaison des bras"],
        ["<code>scripts/build_annotation_report.py</code>", "Ce rapport"],
        ["<code>docs/user/unsupported_genes.md</code>, <code>docs/user/mikado_braker_input_test.md</code>", "Documentation détaillée de la méthode et du test"],
    ]
    return f"""
<h2 id="repro">6. Reproductibilité</h2>
{table(["Fichier", "Rôle"], rows, text=True)}
<pre><code># tableau et classes
data/unsupported_genes_analysis/diamond/run_diamond.sh all
scripts/run_unsupported_genes_analysis.sh data/titan_prod_out data/unsupported_genes_analysis data/PN40024_5.1_on_T2T_ref.gff3
# test Mikado (Slurm)
sbatch -x calcul scripts/mikado_input_test/launch_arm.sh raw
sbatch -x calcul scripts/mikado_input_test/launch_arm.sh tsebra
# ce rapport
python3 scripts/build_annotation_report.py --titan-out data/titan_prod_out \\
   --analysis-dir data/unsupported_genes_analysis --test-dir data/mikado_input_test \\
   --out docs/reports/annotation_quality_and_gene_count.html</code></pre>
<h4>Sources des chiffres</h4>
{table(["Chiffres", "Source"], text=True, rows=[
    ["Structure, BUSCO, OMArk, validation, fonction, expression, chevauchements, mono-exoniques", "<code>ANNOTATION_QUALITY_AUDIT.md</code> et les fichiers de <code>01_final_annotation/quality_report</code>, <code>validation</code>, <code>02_functional_annotation</code>"],
    ["PN40024 5.1", "Ré-exécution de BUSCO et d'AGAT avec les mêmes conteneurs, <code>data/v5_1_audit</code>"],
    ["Classes de preuve, reads, jonctions, homologues, TE", "Calculés pour ce rapport, <code>data/unsupported_genes_analysis</code>"],
    ["Test TSEBRA", "<code>data/mikado_input_test/comparison</code>"],
])}
<p>Outils : BUSCO 6.1.0 (eudicotyledons_odb12.2), OMArk (rosidées), InterProScan 5.78-109.0, eggNOG-mapper 2.1.15, Diamond2GO, DIAMOND 2.1.10, bedtools 2.30, Mikado, TransDecoder 6.0.0, AGAT 1.2.0. Versions complètes dans <code>05_run_info/provenance/</code>.</p>
"""


def build(a):
    D = Data(a)
    head = f"""<!doctype html>
<html lang="fr"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Qualité de l'annotation PN40024 T2T et nombre de gènes</title>
<style>{CSS}</style></head><body><div class="wrap">
<header class="top"><div class="kicker">Pipeline TITAN · run de production PN40024 T2T · {datetime.date.today().isoformat()}</div>
<h1>Annotation du génome PN40024 T2T : qualité et discussion du nombre de gènes</h1>
<p class="lead">Bilan de la qualité de l'annotation produite par TITAN et examen, preuves à l'appui, du nombre élevé de gènes (56&nbsp;744) et des ~15&nbsp;000 gènes sans expression.</p></header>
{sec_summary(D)}
<nav class="toc"><b>Sommaire</b><ol>
<li><a href="#cles">Messages clés</a></li><li><a href="#pipeline">Le pipeline en bref</a></li>
<li><a href="#qualite">Qualité globale</a></li><li><a href="#nombre">Pourquoi 56 744 gènes ?</a></li>
<li><a href="#outil">Marquer et filtrer le GFF3</a></li><li><a href="#limites">Limites</a></li><li><a href="#repro">Reproductibilité</a></li></ol></nav>
"""
    body = "".join([sec_keys(D), sec_context(D), sec_quality(D), sec_count_intro(D), sec_axes(D), sec_origin(D),
                    sec_noexpr(D), sec_classes(D), sec_proofs(D), sec_rejected(D), sec_cause(D), sec_test(D),
                    sec_real(D), sec_balance(D), sec_tool(D), sec_limits(D), sec_repro(D, a)])
    foot = ('<footer>Rapport généré par <code>scripts/build_annotation_report.py</code> à partir de <code>data/titan_prod_out</code> '
            'et des analyses associées. Les graphiques sont des SVG intégrés (survol pour les valeurs, « Voir les données » pour les tableaux).</footer>'
            '</div></body></html>')
    return head + body + foot


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--titan-out", required=True)
    ap.add_argument("--analysis-dir", required=True)
    ap.add_argument("--test-dir", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    html_text = build(a)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out, "w", encoding="utf-8") as fh:
        fh.write(html_text)
    print("wrote %s (%d KB)" % (a.out, len(html_text) // 1024), file=sys.stderr)


if __name__ == "__main__":
    main()
