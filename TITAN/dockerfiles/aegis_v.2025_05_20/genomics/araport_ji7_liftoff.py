#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Jul 03 2023

Program that calls liftoff to generate a range of output
gff files for PN40024 v5 

@author: David Navarro
"""

from modules.tools import bash_run
import os

"""Main Program"""
# INPUT

input_gff_folder = "../../genomes_and_annotation/thale_cress/TAIR10_annotation/"
input_genome_folder = "../../genomes_and_annotation/thale_cress/TAIR10_genome/"
output_gff_folder = "../../genomes_and_annotation/snapdragon/ji7_assembly_3/ji7_annotation_transfer/"
output_genome_folder = "../../genomes_and_annotation/snapdragon/ji7_assembly_3/ji7_genome/"

feature_folder = output_gff_folder + "liftoff_feature_types/"
unmapped_folder = output_gff_folder + "unmapped/"
os.system(f"mkdir -p {unmapped_folder}")

input_gff_files = ["Araport11_annotation_tidy.gff3"]

original_genome_files = ["TAIR10_genome.fasta"]

input_gff_tags = ["Araport11"]

target_tags = ["ji7_3"]
target_genome_fasta_files = ["ama_ji7_genome.fasta"]

"""
Exploring liftoff options:
    
typical command:

liftoff -g input_gff3 -o output_gff3 
        -f types_file target_genome_fasta original_genome_fasta
        
but there are extra options
"""

def_options = ["", "-copies"]

options = ["-mismatch M=1 -gap_open GO=1 -gap_extend GE=0.5", "-copies -mismatch M=1 -gap_open GO=1 -gap_extend GE=0.5"]


liftoff_tags = ["default_liftoff", "copies_liftoff"]
ext = ".gff3"

for n, file in enumerate(input_gff_files):
    print(f"Processing liftoffs for {input_gff_tags[n]}")
    for o, opt in enumerate(options):
        output_processes = []
        for t, target in enumerate(target_genome_fasta_files):
            unmapped_file = (f"{unmapped_folder}{input_gff_tags[n]}_"
                             f"{liftoff_tags[o]}_{target_tags[t]}"
                             "_unmapped.txt")
            print(f"Running a {liftoff_tags[o]} of {input_gff_tags[n]} against"
                  f" {target_tags[t]} target genome\n")
            
            command = (f"liftoff {opt} -g {input_gff_folder}{file} -u {unmapped_file} "
                       f"-dir {output_gff_folder}intermediate_files "
                       f"-f {feature_folder}{input_gff_tags[n]}_types.txt "
                       f"{output_genome_folder}{target} {input_genome_folder}{original_genome_files[n]}")
        

            output_processes.append(bash_run(command, internal_output=True,
                                             tag = (f"{input_gff_tags[n]}_"
                                                    f"{liftoff_tags[o]}_"
                                                    f"{target_tags[t]}_raw"),
                                             folder_out=output_gff_folder,
                                             extension=ext,
                                             standard_output=True,
                                             error_output=True))