#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Processing hop genome and annotation (dovetail) files into some reasonable 
format.

David Navarro Payá
Version 0.0
Date: 01/12/2023
"""

unmasked_genome_full_f = (
    "../../genomes_and_annotation/hop/dovetail_cascade_genome/unmasked/"
    "dovetailCascadeFullAssemblyUnmasked.fasta"
)

unmasked_genome_out_f = (
    "../../genomes_and_annotation/hop/dovetail_cascade_genome/"
    "hop_unmasked_genome.fasta"
)

maker_f = (
    "../../genomes_and_annotation/hop/dovetail_cascade_annotation/"
    "hop/original/makerGenes.gff3"
)

transdecoder_f = (
     "../../genomes_and_annotation/hop/dovetail_cascade_annotation/"
    "hop/original/"
    "transcripts.fasta.transdecoder.genomeCentric.gff3"
)

ids_f = (
     "../../genomes_and_annotation/hop/dovetail_cascade_annotation/"
    "hop/original/combinedGeneModels.txt"
)

gff_f = (
     "../../genomes_and_annotation/hop/dovetail_cascade_annotation/"
    "hop/hop_annotation.gff3"
)

input_gffs = [maker_f, transdecoder_f]

scaffold = ""
seq = ""
scaffolds = {}
# measuring chromosome lengths
f_in = open(unmasked_genome_full_f, "r", encoding="utf-8")
for line in f_in:
    if line.startswith(">"):
        if scaffold != "":
            scaffolds[scaffold] = len(seq)
        scaffold = line[1:-1]
        seq = ""
    else:
        seq += line.strip()
f_in.close()

print(f"There are {len(scaffolds)} scaffolds in the unmasked genome")

# Not making a chromosome 00 just yet since that would be complicated (it would
# change all gff coordinates)

# chromosome sizes in mb:
chr_sizes = {
    476: "chr01",
    434: "chr02",
    370: "chr03",
    345: "chr04",
    328: "chr05",
    317: "chr06",
    304: "chr07",
    291: "chr08",
    185: "chr09",
    424: "chrX",
}

scaffold_to_chr = {}

for size, chr in chr_sizes.items():
    for key, value in scaffolds.items():
        rounded = round(value / 1000000)
        if rounded == size:
            scaffold_to_chr[key] = chr_sizes[rounded]

print(scaffold_to_chr)

f_out = open(unmasked_genome_out_f, "w", encoding="utf-8")
f_out.write("")
f_out.close()

f_out = open(unmasked_genome_out_f, "a", encoding="utf-8")

for scf in scaffold_to_chr:
    f_in = open(unmasked_genome_full_f, "r", encoding="utf-8")
    for line in f_in:
        line = line.strip()
        if line.startswith(">"):
            scaffold = line[1:]
            if scaffold == scf:
                f_out.write(f">{scaffold_to_chr[scf]}\n")
        elif scaffold == scf:
            f_out.write(f"{line}\n")
    f_in.close()

f_in = open(unmasked_genome_full_f, "r", encoding="utf-8")
for line in f_in:
    line = line.strip()
    if line.startswith(">"):
        scaffold = line[1:]
        if scaffold not in scaffold_to_chr:
            f_out.write(f">{scaffold}\n")
    elif scaffold not in scaffold_to_chr:
        f_out.write(f"{line}\n")
f_in.close()
f_out.close()

# reading ID conversions
ids_d = {}

f_in = open(ids_f, "r", encoding="utf-8")
for line in f_in:
    line = line.strip().split("\t")
    unified = line[2].split(".")[0]
    old = line[1].split(".")
    if len(old) == 2:
        old = old[0]
    else:
        old = f"{old[0]}.{old[1]}"

    if old not in ids_d:
        ids_d[old] = unified
f_in.close()

print(f"There are {len(ids_d)} genes in total")

features = ["CDS", "exon", "five_prime_UTR", "gene", "mRNA", "three_prime_UTR"]
sources = ["transdecoder", "maker"]

f_out = open(gff_f, "w", encoding="utf-8")
f_out.write("##gff-version 3\n")
f_out.close()

f_out = open(gff_f, "a", encoding="utf-8")

# parent and ID attributes are updated to their unified IDs
# the way of referring to transcripts should also be fixed at the same time

for f in input_gffs:
    f_in = open(f, "r", encoding="utf-8")
    for line in f_in:
        line = line.strip().split("\t")
        if len(line) > 1:
            source = line[1]
            ft = line[2]
            if ft in features and source in sources:
                scf = line[0]
                att = line[-1].split(";")
                parents_final = []
                ID = ""
                if source == "maker":
                    for a in att:
                        if "ID=" in a:
                            ID = a.split("=")[1].split(".")
                            rest = ".".join(ID[1:])
                            ID = f"ID={ids_d[ID[0]]}"
                            if rest:
                                ID += f".{rest}"
                        elif "Parent=" in a:
                            parents = a.split("=")[1].split(",")
                            for parent in parents:
                                parent = parent.split(".")
                                rest = ".".join(parent[1:])
                                parent = ids_d[parent[0]]
                                if rest:
                                    parent += f".{rest}"
                                parents_final.append(parent)

                else:
                    for a in att:
                        if "ID=" in a:
                            ID = a.split("=")[1].split(".")
                            if ft != "CDS":
                                rest = ".".join(ID[2:])
                                ID = f"{ID[0]}.{ID[1]}"
                                ID = f"ID={ids_d[ID]}"
                                if rest:
                                    ID += f".{rest}"
                            else:
                                rest = ".".join(ID[3:])
                                ID = f"{ID[1]}.{ID[2]}"
                                ID = f"ID={ids_d[ID]}"
                                if rest:
                                    ID += f".{rest}"
                                ID += ".CDS"
                        elif "Parent=" in a:
                            parents = a.split("=")[1].split(",")
                            for parent in parents:
                                parent = parent.split(".")
                                rest = ".".join(parent[2:])
                                parent = f"{parent[0]}.{parent[1]}"
                                parent = ids_d[parent]
                                if rest:
                                    parent += f".{rest}"
                                parents_final.append(parent)
                if scf in scaffold_to_chr:
                    line[0] = scaffold_to_chr[scf]
                if ID:
                    att[0] = ID
                parents_final = ",".join(parents_final)
                if parents_final:
                    att[1] = f"Parent={parents_final}"
                line[-1] = ";".join(att)
                line = "\t".join(line)
                f_out.write(f"{line}\n")
    f_in.close()

f_out.close()
