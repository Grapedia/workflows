#!/usr/bin/env python3
import argparse
import csv
import sys
from functools import lru_cache
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Carry over old liftoff gene IDs onto reciprocal-best-matching "
        "genes in the final annotation, based on an `aegis overlap` correspondence table."
    )
    parser.add_argument("--overlap-csv", required=True)
    parser.add_argument("--gff3", required=True)
    parser.add_argument("--proteins-all", required=True)
    parser.add_argument("--proteins-main", required=True)
    parser.add_argument("--new-origin", default="final")
    parser.add_argument("--old-origin", default="liftoff")
    parser.add_argument("-o", "--output-gff3", required=True)
    parser.add_argument("--output-proteins-all", required=True)
    parser.add_argument("--output-proteins-main", required=True)
    parser.add_argument("--audit-tsv", required=True)
    return parser.parse_args()


def reciprocal_best_matches(overlap_csv, new_origin, old_origin):
    rows = []
    with open(overlap_csv, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("gene_id_A_origin") == new_origin and row.get("gene_id_B_origin") == old_origin:
                rows.append(row)

    best_for_a, best_for_b = {}, {}
    row_by_pair = {}
    for row in rows:
        a, b = row["gene_id_A"], row["gene_id_B"]
        score = float(row["overlap_score"])
        row_by_pair[(a, b)] = row
        if a not in best_for_a or score > best_for_a[a][1]:
            best_for_a[a] = (b, score)
        if b not in best_for_b or score > best_for_b[b][1]:
            best_for_b[b] = (a, score)

    id_map, match_info = {}, {}
    for a, (b, score) in best_for_a.items():
        b_best = best_for_b.get(b)
        if b_best and b_best[0] == a:
            id_map[a] = b
            match_info[a] = row_by_pair[(a, b)]
    return id_map, match_info


def make_remapper(id_map):
    # A linear scan over id_map (tens of thousands of entries) per token,
    # times every ID/Parent token in a genome-scale GFF3/FASTA (hundreds of
    # thousands), is O(n*m) and takes tens of minutes. The gene IDs Mikado/
    # AEGIS/liftoff produce here never contain "_", so old_id + "_" is a
    # prefix of token exactly when old_id equals the token's segment before
    # its first "_" - look that segment up directly (O(1)) instead of
    # scanning. Falls back to the original linear scan (memoized) if that
    # assumption is ever violated, so results stay correct either way; only
    # the common case is fast.
    has_underscore_keys = any("_" in old_id for old_id in id_map)

    @lru_cache(maxsize=None)
    def linear_scan(token):
        for old_id, new_id in id_map.items():
            if token == old_id or token.startswith(old_id + "_"):
                return new_id + token[len(old_id):]
        return token

    def remap_token(token):
        new_id = id_map.get(token)
        if new_id is not None:
            return new_id
        prefix, sep, suffix = token.partition("_")
        if sep:
            new_prefix = id_map.get(prefix)
            if new_prefix is not None:
                return new_prefix + sep + suffix
        return linear_scan(token) if has_underscore_keys else token

    return remap_token


def rewrite_gff3(gff3_path, output_path, remap_token):
    out_lines = []
    original_gene_ids = []
    for line in Path(gff3_path).read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.startswith("#"):
            out_lines.append(line)
            continue
        fields = line.split("\t")
        if len(fields) != 9:
            out_lines.append(line)
            continue

        items = []
        for item in fields[8].split(";"):
            item = item.strip()
            if not item:
                continue
            if "=" in item:
                key, value = item.split("=", 1)
                items.append([key, value])
            else:
                items.append([None, item])

        changed = False
        for pair in items:
            key, value = pair
            if key == "ID":
                if fields[2] == "gene":
                    original_gene_ids.append(value)
                new_value = remap_token(value)
                if new_value != value:
                    pair[1] = new_value
                    changed = True
            elif key == "Parent":
                parts = value.split(",")
                new_parts = [remap_token(part) for part in parts]
                if new_parts != parts:
                    pair[1] = ",".join(new_parts)
                    changed = True

        if changed:
            fields[8] = ";".join(f"{k}={v}" if k else v for k, v in items)
            out_lines.append("\t".join(fields))
        else:
            out_lines.append(line)

    Path(output_path).write_text("\n".join(out_lines) + "\n", encoding="utf-8")
    return original_gene_ids


def rewrite_fasta(in_path, out_path, remap_token):
    out_lines = []
    for line in Path(in_path).read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            header = line[1:]
            parts = header.split(" ", 1)
            new_id = remap_token(parts[0])
            rest = f" {parts[1]}" if len(parts) > 1 else ""
            out_lines.append(f">{new_id}{rest}")
        else:
            out_lines.append(line)
    Path(out_path).write_text("\n".join(out_lines) + "\n", encoding="utf-8")


def write_audit(audit_path, original_gene_ids, id_map, match_info):
    with open(audit_path, "w", encoding="utf-8") as handle:
        handle.write("original_gene_id\tfinal_gene_id\tdecision\toverlap_score\tmin_gene_percent\tmin_cds_percent\n")
        for gene_id in sorted(set(original_gene_ids)):
            if gene_id in id_map:
                info = match_info[gene_id]
                handle.write(
                    f"{gene_id}\t{id_map[gene_id]}\tcarried_over_from_liftoff\t"
                    f"{info.get('overlap_score', '')}\t{info.get('min_gene_percent', '')}\t{info.get('min_CDS_percent', '')}\n"
                )
            else:
                handle.write(f"{gene_id}\t{gene_id}\tnew_vitvi_id\t\t\t\n")


def main():
    args = parse_args()
    id_map, match_info = reciprocal_best_matches(args.overlap_csv, args.new_origin, args.old_origin)
    remap_token = make_remapper(id_map)

    original_gene_ids = rewrite_gff3(args.gff3, args.output_gff3, remap_token)
    rewrite_fasta(args.proteins_all, args.output_proteins_all, remap_token)
    rewrite_fasta(args.proteins_main, args.output_proteins_main, remap_token)
    write_audit(args.audit_tsv, original_gene_ids, id_map, match_info)

    print(
        f"Carried over {len(id_map)} of {len(set(original_gene_ids))} gene IDs from liftoff "
        f"(reciprocal best matches, overlap_score threshold from `aegis overlap` defaults)",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
