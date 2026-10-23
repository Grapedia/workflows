#!/usr/bin/env python3
"""Write a GFF3 subset whose gene and mRNA lines carry the eggNOG-mapper and InterProScan annotation.

NOT part of the pipeline: an analysis tool (e.g. genes absent from the previous annotation, for Apollo).
Genes are kept when listed in --genes (one ID per line), or, with --liftoff-correspondence, when their
decision is new_vitvi_id (no confident Liftoff counterpart in the previous annotation). Both tables are keyed
on gene IDs (TITAN functional outputs of the main proteins). Added attributes:

  description           eggNOG Description, else the InterPro entry names, else the PANTHER / Pfam
                        signature names, else "uncharacterized protein"
  description_source    eggNOG | InterProScan | none
  eggnog_name, eggnog_description, eggnog_cog, eggnog_kegg_ko, eggnog_ec, eggnog_pfams, eggnog_go
  interpro              InterPro entries (IPRxxxxxx:name)
  Dbxref                InterPro:IPR..., Pfam:PF..., PANTHER:PTHR...
  Ontology_term         GO terms from InterProScan
  interpro_signatures   informative signatures (database:accession:name); MobiDBLite, Coils and
                        AntiFam are left out

Values are percent-encoded for GFF3 (; = & , and %). Also writes <out>.tsv with one row per gene.
"""
import argparse
import csv
import re
from collections import OrderedDict

UNINFORMATIVE = {"MobiDBLite", "Coils", "AntiFam"}


def enc(value):
    return (str(value).replace("%", "%25").replace(";", "%3B").replace("=", "%3D")
            .replace("&", "%26").replace(",", "%2C").replace("\t", " "))


def clean(value):
    return "" if value in ("", "-") else value


def load_eggnog(path):
    out = {}
    with open(path) as fh:
        header = None
        for line in fh:
            if line.startswith("##"):
                continue
            f = line.rstrip("\n").split("\t")
            if line.startswith("#query"):
                header = [h.lstrip("#") for h in f]
                continue
            if header:
                out[f[0]] = dict(zip(header, f))
    return out


def load_interproscan(path):
    out = {}
    with open(path) as fh:
        for f in csv.reader(fh, delimiter="\t"):
            if len(f) < 13:
                continue
            gene, analysis, acc, sig_desc = f[0], f[3], f[4], f[5]
            ipr, ipr_desc = f[11], f[12]
            go = f[13] if len(f) > 13 else "-"
            d = out.setdefault(gene, {"ipr": OrderedDict(), "sigs": OrderedDict(), "go": OrderedDict(),
                                      "panther": [], "pfam": []})
            if clean(ipr):
                d["ipr"].setdefault(ipr, clean(ipr_desc))
            if analysis not in UNINFORMATIVE:
                d["sigs"].setdefault(f"{analysis}:{acc}", clean(sig_desc))
                if analysis == "PANTHER" and clean(sig_desc):
                    d["panther"].append(sig_desc)
                if analysis == "Pfam" and clean(sig_desc):
                    d["pfam"].append(sig_desc)
            for term in clean(go).split("|") if clean(go) else []:
                d["go"].setdefault(re.sub(r"\(.*\)$", "", term), None)
    return out


def attributes(gene, eggnog, ips):
    e, i = eggnog.get(gene, {}), ips.get(gene)
    attrs = OrderedDict()
    if clean(e.get("Description", "")):
        attrs["description"], source = e["Description"], "eggNOG"
    elif i and any(i["ipr"].values()):
        attrs["description"], source = " / ".join(v for v in i["ipr"].values() if v), "InterProScan"
    elif i and (i["panther"] or i["pfam"]):
        attrs["description"] = " / ".join(OrderedDict.fromkeys(i["panther"] + i["pfam"]))
        source = "InterProScan"
    else:
        attrs["description"], source = "uncharacterized protein", "none"
    attrs["description_source"] = source
    for key, col in (("eggnog_name", "Preferred_name"), ("eggnog_description", "Description"),
                     ("eggnog_cog", "COG_category"), ("eggnog_kegg_ko", "KEGG_ko"), ("eggnog_ec", "EC"),
                     ("eggnog_pfams", "PFAMs"), ("eggnog_go", "GOs")):
        if clean(e.get(col, "")):
            attrs[key] = e[col].split(",")
    if i:
        if i["ipr"]:
            attrs["interpro"] = [f"{k}:{v}" if v else k for k, v in i["ipr"].items()]
        dbx = [f"InterPro:{k}" for k in i["ipr"]]
        dbx += [f"{s.split(':', 1)[0]}:{s.split(':', 1)[1]}" for s in i["sigs"] if s.split(":")[0] in ("Pfam", "PANTHER")]
        if dbx:
            attrs["Dbxref"] = list(OrderedDict.fromkeys(dbx))
        if i["go"]:
            attrs["Ontology_term"] = list(i["go"])
        if i["sigs"]:
            attrs["interpro_signatures"] = [f"{k}:{v}" if v else k for k, v in i["sigs"].items()]
    return attrs


def render(attrs):
    parts = []
    for k, v in attrs.items():
        parts.append(f"{k}={','.join(enc(x) for x in v)}" if isinstance(v, list) else f"{k}={enc(v)}")
    return ";".join(parts)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gff3", required=True)
    ap.add_argument("--eggnog", required=True, help="emapper.annotations of the main proteins")
    ap.add_argument("--interproscan", required=True, help="InterProScan TSV of the main proteins")
    sel = ap.add_mutually_exclusive_group(required=True)
    sel.add_argument("--genes", help="file with the gene IDs to keep")
    sel.add_argument("--liftoff-correspondence", help="liftoff_gene_id_correspondence.tsv: keep new_vitvi_id genes")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    if args.genes:
        keep = {l.strip() for l in open(args.genes) if l.strip()}
    else:
        with open(args.liftoff_correspondence) as fh:
            keep = {r["final_gene_id"] for r in csv.DictReader(fh, delimiter="\t") if r["decision"] == "new_vitvi_id"}
    eggnog, ips = load_eggnog(args.eggnog), load_interproscan(args.interproscan)

    gene_attrs, n_genes, current = {}, 0, None
    with open(args.gff3) as fh, open(args.out, "w") as out, open(re.sub(r"\.gff3?$", "", args.out) + ".tsv", "w") as tsv:
        out.write("##gff-version 3\n")
        tsv.write("gene_id\tchrom\tstart\tend\tstrand\tdescription\tdescription_source\teggnog_name\t"
                  "eggnog_kegg_ko\tinterpro\n")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) != 9:
                continue
            if f[2] == "gene":
                gid = re.search(r"ID=([^;]+)", f[8]).group(1)
                current = gid if gid in keep else None
                if current:
                    n_genes += 1
                    gene_attrs[gid] = attributes(gid, eggnog, ips)
                    a = gene_attrs[gid]
                    if n_genes > 1:
                        out.write("###\n")
                    tsv.write("\t".join([gid, f[0], f[3], f[4], f[6], a["description"], a["description_source"],
                                         ",".join(a.get("eggnog_name", [])), ",".join(a.get("eggnog_kegg_ko", [])),
                                         " | ".join(a.get("interpro", []))]) + "\n")
            if current is None:
                continue
            if f[2] in ("gene", "mRNA"):
                f[8] = f[8].rstrip(";") + ";" + render(gene_attrs[current])
            out.write("\t".join(f) + "\n")
        if n_genes:
            out.write("###\n")
    print(f"genes written: {n_genes} (requested: {len(keep)})")


if __name__ == "__main__":
    main()
