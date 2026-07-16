#!/usr/bin/env python3
import argparse
import csv
import gzip
import re
import sys
from pathlib import Path


SAMPLE_ID_RE = re.compile(r"^[A-Za-z0-9_.-]+$")
RNA_COLUMNS = ["sampleID", "SRA_or_FASTQ", "library_layout"]
PROTEIN_COLUMNS = ["organism", "filename"]
RNA_SOURCES = {"SRA", "FASTQ", "FASTA"}
RNA_LAYOUTS = {"single", "paired", "long"}
EGAPX_EXECUTORS = {"docker", "singularity", "apptainer"}


def error(message):
    return f"ERROR: {message}"


def read_fasta(path):
    records = {}
    seen_ids = set()
    current = None
    chunks = []
    with path.open() as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current:
                    records[current] = "".join(chunks)
                current = line[1:].split()[0]
                if not current:
                    raise ValueError(f"{path}: line {line_number}: FASTA header has no sequence id")
                if current in seen_ids:
                    raise ValueError(f"{path}: duplicate FASTA sequence id: {current}")
                seen_ids.add(current)
                chunks = []
            else:
                if current is None:
                    raise ValueError(f"{path}: line {line_number}: sequence appears before first FASTA header")
                chunks.append(line)
    if current:
        records[current] = "".join(chunks)
    if not records:
        raise ValueError(f"{path}: FASTA contains no sequences")
    empty = [name for name, seq in records.items() if not seq]
    if empty:
        raise ValueError(f"{path}: empty FASTA sequence(s): {', '.join(empty)}")
    return records


def validate_gff3(path, seqs):
    ids = set()
    parents = []
    feature_count = 0
    with path.open() as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            feature_count += 1
            fields = line.split("\t")
            if len(fields) != 9:
                raise ValueError(f"{path}: line {line_number}: expected 9 GFF3 columns, found {len(fields)}")
            seqid, _, _, start, end, _, _, _, attributes = fields
            if seqid not in seqs:
                raise ValueError(f"{path}: line {line_number}: seqid '{seqid}' is absent from FASTA")
            try:
                start_i = int(start)
                end_i = int(end)
            except ValueError:
                raise ValueError(f"{path}: line {line_number}: start/end must be integers")
            if start_i < 1 or end_i < start_i or end_i > len(seqs[seqid]):
                raise ValueError(f"{path}: line {line_number}: invalid coordinates {start}-{end} for {seqid}")
            attrs = {}
            for item in attributes.split(";"):
                if "=" in item:
                    key, value = item.split("=", 1)
                    attrs[key] = value
            if "ID" in attrs:
                if attrs["ID"] in ids:
                    raise ValueError(f"{path}: line {line_number}: duplicate GFF3 ID '{attrs['ID']}'")
                ids.add(attrs["ID"])
            if "Parent" in attrs:
                parents.extend(parent for parent in attrs["Parent"].split(",") if parent)
    if feature_count == 0:
        raise ValueError(f"{path}: GFF3 contains no features")
    missing = sorted({parent for parent in parents if parent not in ids})
    if missing:
        raise ValueError(f"{path}: missing Parent target(s): {', '.join(missing)}")


def validate_fastq_gz(path):
    with gzip.open(path, "rt") as handle:
        lines = [line.rstrip("\n") for _, line in zip(range(400), handle)]
    if not lines or len(lines) % 4 != 0:
        raise ValueError(f"{path}: FASTQ does not contain complete records")
    for offset in range(0, len(lines), 4):
        name, sequence, plus, quality = lines[offset : offset + 4]
        if not name.startswith("@") or plus != "+" or len(sequence) != len(quality):
            raise ValueError(f"{path}: invalid FASTQ record starting at sampled line {offset + 1}")


def validate_existing_file(path, label, errors):
    if not path.exists() or not path.is_file():
        errors.append(error(f"{label} does not exist or is not a file: {path}"))


def validate_rnaseq_samplesheet(path, rnaseq_data_dir):
    errors = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            return [error(f"{path}: RNA-seq samplesheet is empty")]
        missing = [column for column in RNA_COLUMNS if column not in reader.fieldnames]
        if missing:
            return [error(f"{path}: missing RNA-seq column(s): {', '.join(missing)}")]
        rows = list(reader)

    if not rows:
        errors.append(error(f"{path}: RNA-seq samplesheet contains no data rows"))

    seen = set()
    layouts = set()
    for index, row in enumerate(rows, 2):
        sample_id = (row.get("sampleID") or "").strip()
        source = (row.get("SRA_or_FASTQ") or "").strip().upper()
        layout = (row.get("library_layout") or "").strip().lower()
        if not sample_id:
            errors.append(error(f"{path}: line {index}: sampleID is required"))
            continue
        if not SAMPLE_ID_RE.match(sample_id):
            errors.append(error(f"{path}: line {index}: sampleID '{sample_id}' contains unsupported characters"))
        if sample_id in seen:
            errors.append(error(f"{path}: line {index}: duplicate sampleID '{sample_id}'"))
        seen.add(sample_id)
        if source not in RNA_SOURCES:
            errors.append(error(f"{path}: line {index}: SRA_or_FASTQ must be one of {sorted(RNA_SOURCES)}, got '{source}'"))
        if layout not in RNA_LAYOUTS:
            errors.append(error(f"{path}: line {index}: library_layout must be one of {sorted(RNA_LAYOUTS)}, got '{layout}'"))
        else:
            layouts.add(layout)
        if source == "FASTA" and layout != "long":
            errors.append(error(f"{path}: line {index}: FASTA reads are only supported for library_layout=long"))
        if source == "FASTQ":
            if layout == "paired":
                files = [rnaseq_data_dir / f"{sample_id}_1.fastq.gz", rnaseq_data_dir / f"{sample_id}_2.fastq.gz"]
            elif layout in {"single", "long"}:
                files = [rnaseq_data_dir / f"{sample_id}.fastq.gz"]
            else:
                files = []
            for fastq in files:
                validate_existing_file(fastq, f"RNA-seq FASTQ for sample {sample_id}", errors)
                if fastq.exists():
                    try:
                        validate_fastq_gz(fastq)
                    except Exception as exc:
                        errors.append(error(str(exc)))
        if source == "FASTA" and layout == "long":
            fasta = rnaseq_data_dir / f"{sample_id}.fasta"
            validate_existing_file(fasta, f"RNA-seq FASTA for sample {sample_id}", errors)
            if fasta.exists():
                try:
                    read_fasta(fasta)
                except Exception as exc:
                    errors.append(error(str(exc)))

    if "single" not in layouts and "paired" not in layouts:
        errors.append(error(f"{path}: at least one short-read row with library_layout single or paired is required"))
    return errors


def validate_protein_samplesheet(path, root):
    errors = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            return [error(f"{path}: protein samplesheet is empty")]
        missing = [column for column in PROTEIN_COLUMNS if column not in reader.fieldnames]
        if missing:
            return [error(f"{path}: missing protein column(s): {', '.join(missing)}")]
        rows = list(reader)

    if not rows:
        errors.append(error(f"{path}: protein samplesheet contains no data rows"))

    seen = set()
    for index, row in enumerate(rows, 2):
        organism = (row.get("organism") or "").strip()
        filename = (row.get("filename") or "").strip()
        if not organism:
            errors.append(error(f"{path}: line {index}: organism is required"))
        if not filename:
            errors.append(error(f"{path}: line {index}: filename is required"))
            continue
        protein_path = Path(filename)
        if not protein_path.is_absolute():
            protein_path = root / protein_path
        if filename in seen:
            errors.append(error(f"{path}: line {index}: duplicate protein filename '{filename}'"))
        seen.add(filename)
        validate_existing_file(protein_path, f"protein FASTA for {organism or filename}", errors)
        if protein_path.exists():
            try:
                read_fasta(protein_path)
            except Exception as exc:
                errors.append(error(str(exc)))
    return errors


def parse_simple_yaml(path):
    values = {}
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#") or ":" not in line:
                continue
            key, value = line.split(":", 1)
            values[key.strip()] = value.strip().strip("'\"")
    return values


def validate_egapx_yaml(path):
    errors = []
    values = parse_simple_yaml(path)
    for key in ("taxid", "organism", "genome"):
        if not values.get(key):
            errors.append(error(f"{path}: EGAPx YAML must define '{key}'"))
    genome = values.get("genome")
    if genome:
        genome_path = Path(genome)
        if not genome_path.is_absolute():
            genome_path = path.parent / genome_path
        validate_existing_file(genome_path, "EGAPx genome", errors)
        if genome_path.exists():
            try:
                read_fasta(genome_path)
            except Exception as exc:
                errors.append(error(str(exc)))
    return errors


def validate_args(args):
    errors = []
    required_files = {
        "new assembly": args.new_assembly,
        "previous assembly": args.previous_assembly,
        "previous annotations": args.previous_annotations,
        "RNA-seq samplesheet": args.rnaseq_samplesheet,
        "protein samplesheet": args.protein_samplesheet,
        "EGAPx YAML": args.egapx_paramfile,
    }
    for label, path in required_files.items():
        validate_existing_file(path, label, errors)
    if not args.rnaseq_data_dir.exists() or not args.rnaseq_data_dir.is_dir():
        errors.append(error(f"RNA-seq data directory does not exist or is not a directory: {args.rnaseq_data_dir}"))

    if errors:
        return errors

    try:
        previous_seqs = read_fasta(args.previous_assembly)
    except Exception as exc:
        errors.append(error(str(exc)))
        previous_seqs = {}
    try:
        read_fasta(args.new_assembly)
    except Exception as exc:
        errors.append(error(str(exc)))
    if previous_seqs:
        try:
            validate_gff3(args.previous_annotations, previous_seqs)
        except Exception as exc:
            errors.append(error(str(exc)))

    errors.extend(validate_rnaseq_samplesheet(args.rnaseq_samplesheet, args.rnaseq_data_dir))
    errors.extend(validate_protein_samplesheet(args.protein_samplesheet, args.project_dir))
    errors.extend(validate_egapx_yaml(args.egapx_paramfile))

    if args.egapx_executor not in EGAPX_EXECUTORS:
        errors.append(error(f"egapx_executor must be one of {sorted(EGAPX_EXECUTORS)}, got '{args.egapx_executor}'"))
    for label, value in {"PSICLASS_vd_option": args.psiclass_vd, "PSICLASS_c_option": args.psiclass_c}.items():
        try:
            float(value)
        except ValueError:
            errors.append(error(f"{label} must be numeric, got '{value}'"))
    return errors


def build_parser():
    parser = argparse.ArgumentParser(description="Validate TITAN production inputs before launching Nextflow tasks.")
    parser.add_argument("--project-dir", type=Path, default=Path.cwd())
    parser.add_argument("--new-assembly", type=Path, required=True)
    parser.add_argument("--previous-assembly", type=Path, required=True)
    parser.add_argument("--previous-annotations", type=Path, required=True)
    parser.add_argument("--rnaseq-samplesheet", type=Path, required=True)
    parser.add_argument("--rnaseq-data-dir", type=Path, required=True)
    parser.add_argument("--protein-samplesheet", type=Path, required=True)
    parser.add_argument("--egapx-paramfile", type=Path, required=True)
    parser.add_argument("--egapx-executor", default="docker")
    parser.add_argument("--psiclass-vd", default="5.0")
    parser.add_argument("--psiclass-c", default="0.03")
    return parser


def main():
    args = build_parser().parse_args()
    errors = validate_args(args)
    if errors:
        print("TITAN input validation failed:", file=sys.stderr)
        for item in errors:
            print(f"  {item}", file=sys.stderr)
        return 1
    print("TITAN input validation OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
