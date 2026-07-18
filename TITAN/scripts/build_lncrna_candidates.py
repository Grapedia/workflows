#!/usr/bin/env python3
"""Build preliminary lncRNA candidates from transcript GTF evidence."""

import argparse
import sys
from collections import defaultdict
from pathlib import Path
from urllib.parse import quote


def parse_attributes(raw):
    attributes = {}
    for item in raw.strip().strip(";").split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
        elif " " in item:
            key, value = item.split(" ", 1)
            value = value.strip().strip('"')
        else:
            key, value = item, ""
        attributes.setdefault(key, []).append(value)
    return attributes


def first_attr(attributes, *keys):
    for key in keys:
        values = attributes.get(key, [])
        if values:
            return values[0]
    return ""


def read_fasta(path):
    sequences = {}
    name = None
    parts = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name:
                    sequences[name] = "".join(parts).upper()
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
    if name:
        sequences[name] = "".join(parts).upper()
    return sequences


def intervals_overlap(left, right):
    return left[0] == right[0] and left[1] <= right[2] and right[1] <= left[2]


def load_blocking_intervals(paths, feature_types=None):
    intervals = []
    for path in paths:
        with path.open(encoding="utf-8") as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 9:
                    continue
                if feature_types and fields[2] not in feature_types:
                    continue
                try:
                    intervals.append((fields[0], int(fields[3]), int(fields[4])))
                except ValueError:
                    continue
    return intervals


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1].upper()


def transcript_sequence(genome, seqid, exons, strand):
    sequence = genome.get(seqid, "")
    chunks = []
    for start, end in sorted(exons):
        chunks.append(sequence[start - 1 : end])
    joined = "".join(chunks).upper()
    return reverse_complement(joined) if strand == "-" else joined


def collect_transcripts(gtf_paths):
    transcripts = {}
    exons = defaultdict(list)
    for gtf_path in gtf_paths:
        source_name = gtf_path.name
        with gtf_path.open(encoding="utf-8") as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 9:
                    continue
                seqid, source, feature_type, start, end, score, strand, phase, raw_attrs = fields
                try:
                    start_i, end_i = int(start), int(end)
                except ValueError:
                    continue
                attrs = parse_attributes(raw_attrs)
                if feature_type in {"transcript", "mRNA"}:
                    transcript_id = first_attr(attrs, "transcript_id", "ID")
                    if not transcript_id:
                        continue
                    transcripts.setdefault(
                        transcript_id,
                        {
                            "seqid": seqid,
                            "source": source or source_name,
                            "start": start_i,
                            "end": end_i,
                            "strand": strand if strand in {"+", "-"} else ".",
                        },
                    )
                elif feature_type == "exon":
                    transcript_id = first_attr(attrs, "transcript_id", "Parent")
                    if transcript_id:
                        exons[transcript_id].append((start_i, end_i))
    for transcript_id, record in list(transcripts.items()):
        if not exons[transcript_id]:
            exons[transcript_id].append((record["start"], record["end"]))
        record["exons"] = sorted(exons[transcript_id])
        record["length"] = sum(end - start + 1 for start, end in record["exons"])
    return transcripts


def load_cpat_coding_ids(path, cutoff):
    if not path:
        return set()
    coding_ids = set()
    with path.open(encoding="utf-8") as handle:
        header = []
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            lowered = [field.lower() for field in fields]
            if not header and any("coding" in field or "prob" in field for field in lowered):
                header = lowered
                continue
            probability = None
            if header:
                for key in ("coding_prob", "coding_probability", "coding probability", "prob"):
                    if key in header:
                        try:
                            probability = float(fields[header.index(key)])
                        except (IndexError, ValueError):
                            probability = None
                        break
            if probability is None:
                for value in reversed(fields[1:]):
                    try:
                        probability = float(value)
                        break
                    except ValueError:
                        continue
            if probability is not None and probability >= cutoff:
                coding_ids.add(fields[0])
    return coding_ids


def write_outputs(transcripts, genome, blockers, min_length, prefix, cpat_coding_ids=None):
    gff3_path = Path(f"{prefix}.gff3")
    gtf_path = Path(f"{prefix}.gtf")
    fasta_path = Path(f"{prefix}.fasta")
    summary_path = Path("lncrna_classification_summary.tsv")
    mqc_path = Path("lncrna_candidates_mqc.tsv")

    kept = []
    excluded_short = 0
    excluded_overlap = 0
    excluded_cpat = 0
    for transcript_id, record in sorted(transcripts.items()):
        if record["length"] < min_length:
            excluded_short += 1
            continue
        interval = (record["seqid"], min(start for start, _ in record["exons"]), max(end for _, end in record["exons"]))
        if any(intervals_overlap(interval, blocker) for blocker in blockers):
            excluded_overlap += 1
            continue
        candidate_id = f"lncrna_candidate_{len(kept) + excluded_cpat + 1}.t1"
        if cpat_coding_ids and candidate_id in cpat_coding_ids:
            excluded_cpat += 1
            continue
        kept.append((transcript_id, record, interval))

    with gff3_path.open("w", encoding="utf-8") as gff3, gtf_path.open("w", encoding="utf-8") as gtf, fasta_path.open("w", encoding="utf-8") as fasta:
        gff3.write("##gff-version 3\n")
        for index, (transcript_id, record, interval) in enumerate(kept, 1):
            gene_id = f"lncrna_candidate_{index}"
            encoded_transcript = quote(transcript_id, safe="._:-")
            attrs_gene = f"ID={gene_id};Name={encoded_transcript};biotype=lncRNA_candidate"
            attrs_tx = f"ID={gene_id}.t1;Parent={gene_id};Name={encoded_transcript};biotype=lncRNA_candidate"
            gff3.write("\t".join([record["seqid"], "TITAN_lncRNA", "gene", str(interval[1]), str(interval[2]), ".", record["strand"], ".", attrs_gene]) + "\n")
            gff3.write("\t".join([record["seqid"], "TITAN_lncRNA", "lnc_RNA", str(interval[1]), str(interval[2]), ".", record["strand"], ".", attrs_tx]) + "\n")
            for exon_index, (start, end) in enumerate(record["exons"], 1):
                exon_attrs = f"ID={gene_id}.t1.exon{exon_index};Parent={gene_id}.t1"
                gff3.write("\t".join([record["seqid"], "TITAN_lncRNA", "exon", str(start), str(end), ".", record["strand"], ".", exon_attrs]) + "\n")
                gtf_attrs = f'gene_id "{gene_id}"; transcript_id "{gene_id}.t1"; source_transcript "{transcript_id}";'
                gtf.write("\t".join([record["seqid"], "TITAN_lncRNA", "exon", str(start), str(end), ".", record["strand"], ".", gtf_attrs]) + "\n")
            seq = transcript_sequence(genome, record["seqid"], record["exons"], record["strand"])
            fasta.write(f">{gene_id}.t1 source={encoded_transcript} length={len(seq)}\n")
            for pos in range(0, len(seq), 80):
                fasta.write(seq[pos : pos + 80] + "\n")

    with summary_path.open("w", encoding="utf-8") as summary:
        summary.write("class\tcount\n")
        summary.write(f"lncRNA_candidate\t{len(kept)}\n")
        summary.write(f"excluded_short\t{excluded_short}\n")
        summary.write(f"excluded_overlap_coding_or_ncrna\t{excluded_overlap}\n")
        summary.write(f"excluded_cpat_coding\t{excluded_cpat}\n")

    with mqc_path.open("w", encoding="utf-8") as mqc:
        mqc.write("# id: titan_lncrna_candidates\n")
        mqc.write("# section_name: 'TITAN lncRNA candidates'\n")
        mqc.write("# description: 'Preliminary lncRNA candidates after length and coding/ncRNA overlap filters.'\n")
        mqc.write("# plot_type: 'table'\n")
        mqc.write("Metric\tCount\n")
        mqc.write(f"lncRNA candidates\t{len(kept)}\n")
        mqc.write(f"Excluded short transcripts\t{excluded_short}\n")
        mqc.write(f"Excluded coding/ncRNA overlaps\t{excluded_overlap}\n")
        mqc.write(f"Excluded CPAT-plant coding calls\t{excluded_cpat}\n")

    return len(kept)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genome", required=True, type=Path)
    parser.add_argument("--final-annotation", required=True, type=Path)
    parser.add_argument("--trna-gff3", required=True, type=Path)
    parser.add_argument("--rfam-gff3", required=True, type=Path)
    parser.add_argument("--min-length", required=True, type=int)
    parser.add_argument("--cpat-best-tsv", type=Path)
    parser.add_argument("--cpat-cutoff", type=float, default=0.46)
    parser.add_argument("--output-prefix", default="lncrna_candidates")
    parser.add_argument("gtf", nargs="+", type=Path)
    args = parser.parse_args(argv)

    genome = read_fasta(args.genome)
    transcripts = collect_transcripts(args.gtf)
    blockers = load_blocking_intervals([args.final_annotation], {"CDS"})
    blockers.extend(load_blocking_intervals([args.trna_gff3, args.rfam_gff3]))
    cpat_coding_ids = load_cpat_coding_ids(args.cpat_best_tsv, args.cpat_cutoff)
    write_outputs(transcripts, genome, blockers, args.min_length, args.output_prefix, cpat_coding_ids)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
