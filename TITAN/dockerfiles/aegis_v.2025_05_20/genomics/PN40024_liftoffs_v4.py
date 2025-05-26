#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 13:07:45 2023

Program that calls liftoff to generate a range of output
gff files for PN40024 v4 which has genome divided in ref and alt

@author: David Navarro
"""

from modules.tools import bash_run
import os

"""Main Program"""
# INPUT
root = "../../genomes_and_annotation/grapevine/PN40024/"

output_folder = root + "40X_annotation_transfer/"
feature_folder = output_folder + "liftoff_feature_types/"
unmapped_folder = output_folder + "unmapped/"
os.system(f"mkdir -p {unmapped_folder}")
os.system(f"mkdir -p {output_folder}original/")

input_gff_files = ["8X_annotation_genoscope/8x_tidy.gff3",
                   "12X_annotation_NCBI/NCBI_original_tidy.gff3",
                   "12Xv2_annotation_transfer_URGI/v0_on_12Xv2_tidy.gff3",
                   "12Xv2_annotation_transfer_URGI/v1_on_12Xv2_tidy.gff3",
                   "12Xv2_annotation_transfer_URGI/v2_on_12Xv2_tidy.gff3",
                   "12Xv2_annotation_VCostv3/v3_tidy.gff3"]

original_genome_files = ["8X_genome_genoscope/8x.fasta",
                         "12X_genome_NCBI/12X_genome_NCBI.fasta",
                         "12Xv2_genome/12Xv2_genome.fasta",
                         "12Xv2_genome/12Xv2_genome.fasta",
                         "12Xv2_genome/12Xv2_genome.fasta",
                         "12Xv2_genome/12Xv2_genome.fasta"]

input_gff_tags = ["8x", "NCBI_original", "v0_on_12Xv2", "v1_on_12Xv2", 
                  "v2_on_12Xv2", "v3"]

target_tags = ["ref", "alt"]
target_genome_fasta_files = ["40X_genome/v4_genome_ref.fasta",
                             "40X_genome/v4_genome_alt.fasta"]

"""
Exploring liftoff options:
    
typical command:

liftoff -g input_gff3 -o output_gff3 
        -f types_file target_genome_fasta original_genome_fasta
        
but there are extra options
"""

options = ["", "-copies"]
liftoff_tags = ["default_liftoff", "copies_liftoff"]
ext = ".gff3"

for n, file in enumerate(input_gff_files):
    print(f"Processing liftoffs for {input_gff_tags[n]}")
    for o, opt in enumerate(options):
        output_processes = []
        output_unmapped = []
        liftoff_file_1 = (f"{input_gff_tags[n]}_{liftoff_tags[o]}_"
                          f"{target_tags[0]}_raw{ext}")
        liftoff_file_2 = (f"{input_gff_tags[n]}_{liftoff_tags[o]}_"
                          f"{target_tags[1]}_raw{ext}")
        if (os.path.isfile(f"{output_folder}original/{liftoff_file_1}") and
            os.path.isfile(f"{output_folder}original/{liftoff_file_2}")):
            continue
        for t, target in enumerate(target_genome_fasta_files):
            unmapped_file = (f"{unmapped_folder}{input_gff_tags[n]}_"
                             f"{liftoff_tags[o]}_{target_tags[t]}"
                             "_unmapped.txt")
            print(f"Running a {liftoff_tags[o]} of {input_gff_tags[n]} against"
                  f" {target_tags[t]} target genome\n")
            
            command = (f"liftoff {opt} -g {root}{file} -u {unmapped_file} "
                       f"-dir {output_folder}intermediate_files "
                       f"-f {feature_folder}{input_gff_tags[n]}_types.txt "
                       f"{root}{target} {root}{original_genome_files[n]}")
        

            output_processes.append(bash_run(command, internal_output=True,
                                             tag = (f"{input_gff_tags[n]}_"
                                                    f"{liftoff_tags[o]}_"
                                                    f"{target_tags[t]}_raw"),
                                             folder_out=output_folder,
                                             extension=ext,
                                             standard_output=True,
                                             error_output=True))
            
            f = open(unmapped_file, "r", encoding="utf-8")
            output_unmapped.append(f.read()[:-1].split("\n"))
            f.close()
        
        # Unmapped genes in both ref and alt liftoffs
        unmapped_text = ""
        unmapped_file = (f"{unmapped_folder}{input_gff_tags[n]}_"
                         f"{liftoff_tags[o]}_merged_unmapped.txt")
        for gene in output_unmapped[0]:
            if gene in output_unmapped[1]:
                unmapped_text += gene + "\n"
        
        f = open(unmapped_file, "w", encoding="utf-8")
        f.write(unmapped_text)
        f.close()

for n, file in enumerate(input_gff_files):
    print(f"Fixing liftoff format for  {input_gff_tags[n]}")
    for o, opt in enumerate(options):
        for t, target in enumerate(target_genome_fasta_files):
            f_tag = f"{input_gff_tags[n]}_{liftoff_tags[o]}_{target_tags[t]}"
            out = ""
            try:
                os.system(f"mv {output_folder}{input_gff_tags[n]}_"
                        f"{liftoff_tags[o]}_{target_tags[t]}_raw{ext} "
                        f"{output_folder}original/{input_gff_tags[n]}_"
                        f"{liftoff_tags[o]}_{target_tags[t]}_raw{ext}")
            except:
                print(f"raw {output_folder}{input_gff_tags[n]}_"
                      f"{liftoff_tags[o]}_{target_tags[t]}_raw{ext} file "
                      "already moved")
            f = open(f"{output_folder}original/{f_tag}_raw{ext}", "r", 
                     encoding="utf-8")
            header = []
            for line in f:
                if not line:
                    continue
                if line.startswith("#"):
                    if line == "###\n":
                        out += line
                    elif (line.startswith("##gff")
                          or line.startswith("# Liftoff")):
                        if line not in header:
                            header.append(line)
                    continue
                
                line = line.strip().split("\t")
                if len(line) <= 2:
                    continue
                else:
                    out += "\t".join(line) + "\n"

            f.close()

            out =  "".join(header) + out

            f = open(f"{output_folder}{f_tag}.gff3",
                     "w", encoding="utf-8")
            f.write(out)
            f.close()

# Assuming two target tags, one ref, one alt
for n, file in enumerate(input_gff_files):
    print(f"Merging liftoffs for  {input_gff_tags[n]}")
    for o, opt in enumerate(options):
        f_tag_1 = f"{input_gff_tags[n]}_{liftoff_tags[o]}_{target_tags[0]}"
        f_tag_1 += f"{ext}"
        f_tag_2 = f"{input_gff_tags[n]}_{liftoff_tags[o]}_{target_tags[1]}"
        f_tag_2 += f"{ext}"
        header = []
        out = ""
        f_in = open(f"{output_folder}{f_tag_1}", "r", encoding="utf-8")
        for line in f_in:
            if not line:
                continue
            if line == ("###\n"):
                out += line
                continue
            if line.startswith("#"):
                if line not in header:
                    header.append(line)
            else:
                out += line
        f_in.close()

        f_in = open(f"{output_folder}{f_tag_2}", "r", encoding="utf-8")
        for line in f_in:
            if not line:
                continue
            if line == ("###\n"):
                out += line
                continue
            if line.startswith("#"):
                if line not in header:
                    header.append(line)
            else:
                out += line
        f_in.close()

        out =  "".join(header) + out
        out_tag = f"{input_gff_tags[n]}_{liftoff_tags[o]}_merged.gff3"
        f_out = open(f"{output_folder}{out_tag}", "w", encoding="utf-8")
        f_out.write(out)
        f_out.close()
