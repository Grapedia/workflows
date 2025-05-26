#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on May 8 2025

@author: David Navarro
"""

import os

"""Main Program"""
# INPUT
root = "../../genomes_and_annotation/"

root = f"{os.path.abspath(root)}/"

output_folder = root + "grapevine/PN40024/T2T_annotation_transfer/other_species/"
feature_folder = output_folder + "liftoff_feature_types/"
unmapped_folder = output_folder + "unmapped/"
os.system(f"mkdir -p {unmapped_folder}")
os.system(f"mkdir -p {output_folder}original/")

input_gff_files = ["snapdragon/ji7_assembly_3/aegis_output/gffs/ji7_3_on_ji7_3_aegis_fcounts.gff3"]

original_genome_files = ["snapdragon/ji7_assembly_3/ji7_genome/ama_ji7_genome.fasta"]

input_gff_tags = ["ji7"]

target_tag = "T2T_ref"
target_genome_fasta_file = "grapevine/PN40024/T2T_genome/T2T_ref.fasta"

liftoff = True
lifton = False

"""
Exploring liftoff options:
    
typical command:

liftoff -g input_gff3 -o output_gff3 -f types_file target_genome_fasta original_genome_fasta
        
but there are extra options
"""

# change to empty string if you want to use the native envitonment's lifton installation

singularity_prefix = "singularity run -B /storage/tom:/storage/tom --no-home /storage/tom/lifton.sif"

options = ["", "-copies", "-copies -flank 0.1"]
liftoff_tags = ["default_liftoff", "copies_liftoff", "copies_flank_liftoff"]
lifton_tags = ["default_lifton", "copies_lifton", "copies_flank_lifton"]

cores = 8

# options = [""]
# liftoff_tags = ["default_liftoff"]
# lifton_tags = ["default_lifton"]

ext = ".gff3"

for n, file in enumerate(input_gff_files):
    if lifton and liftoff:
        print(f"Processing liftoffs and liftons for {input_gff_tags[n]}")
    elif lifton:
        print(f"Processing liftons for {input_gff_tags[n]}")
    else:
        print(f"Processing liftoffs for {input_gff_tags[n]}")
    for o, opt in enumerate(options):

        if liftoff:

            liftoff_tag = f"{input_gff_tags[n]}_{liftoff_tags[o]}_{target_tag}"
            liftoff_file = f"{liftoff_tag}_raw{ext}"

            liftoff_file_path = f"{output_folder}original/{liftoff_file}"

            unmapped_file = f"{unmapped_folder}{liftoff_tag}_unmapped.txt"

            intermediate_path = f"{output_folder}intermediate_files"

            if os.path.isdir(intermediate_path):
                os.system(f"rm -r {intermediate_path}")

            if not os.path.isfile(liftoff_file_path):
                print(f"Running a {liftoff_tags[o]} of {input_gff_tags[n]} against {target_tag} target genome\n")
                
                command = (f"liftoff {opt} -g {root}{file} -o {liftoff_file_path} -u {unmapped_file} -dir {intermediate_path} -f {feature_folder}{input_gff_tags[n]}_types.txt {root}{target_genome_fasta_file} {root}{original_genome_files[n]}")
            
                os.system(command)

                if os.path.isdir(intermediate_path):
                    os.system(f"rm -r {intermediate_path}")

            if not os.path.isfile(f"{output_folder}{liftoff_tag}.gff3"):
                f = open(f"{liftoff_file_path}", "r", encoding="utf-8")
                header = []
                out = []
                for line in f:
                    if not line:
                        continue
                    if line.startswith("#"):
                        if line == "###\n":
                            out.append(line)
                        elif line.startswith("##gff") or line.startswith("# Liftoff"):
                            if line not in header:
                                header.append(line)
                        continue
                    
                    temp = line.strip().split("\t")
                    if len(temp) <= 2:
                        continue

                    out.append(line)

                f.close()

                out = header + out

                out = "".join(out)

                f = open(f"{output_folder}{liftoff_tag}.gff3", "w", encoding="utf-8")
                f.write(out)
                f.close()

        if lifton:

            lifton_tag = f"{input_gff_tags[n]}_{lifton_tags[o]}_{target_tag}"

            lifton_file = f"{lifton_tag}{ext}"

            lifton_file_path = f"{output_folder}{lifton_file}"

            unmapped_file = f"{unmapped_folder}{lifton_tag}_unmapped.txt"

            lifton_output_path = f"{output_folder}lifton_output/"

            if os.path.isdir(lifton_output_path):
                os.system(f"rm -r {lifton_output_path}")

            if not os.path.isfile(lifton_file_path):

                print(f"Running a {lifton_tags[o]} of {input_gff_tags[n]} against {target_tag} target genome\n")
                
                command = (f"{singularity_prefix} lifton {opt} -t {cores} -g {root}{file} -o {lifton_file_path} -u {unmapped_file} -f {feature_folder}{input_gff_tags[n]}_types.txt {root}{target_genome_fasta_file} {root}{original_genome_files[n]}")
            
                os.system(command)

                if os.path.isdir(lifton_output_path):
                    os.system(f"rm -r {lifton_output_path}")


