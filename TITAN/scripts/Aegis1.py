#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import time
import os
import argparse
from modules.tools import pickle_load, pickle_save
from modules.annotation import Annotation
from modules.genome import Genome

# Main function
def main(args):
    
    chosen_chromosomes = [args.chromosome] if args.chromosome else None
    if chosen_chromosomes:
        print(f"-------------------- Chosen chromosome is : {chosen_chromosomes}")
    else:
        print("-------------------- Aegis launched on all chromosomes.")

    print("-------------------- Starting script execution")

    # Verify if outdir exists or create it
    if not os.path.exists(args.output_dir):
        print(f"Output directory {args.output_dir} does not exist. Creating it.")
        os.makedirs(args.output_dir)

    # Upload the genome assembly
    print(f"-------------------- Loading genome: {args.genome_name} from {args.genome_path}")
    genome = Genome(args.genome_name, args.genome_path)

    # define evidences
    print("-------------------- Defining evidences")
    evidences = {
        'ab_initio': [
            ('augustus', args.augustus_path),
            ('genemark', args.genemark_path),
            ('liftoff', args.liftoff_path)
        ],
        'transcriptome': [
            ('psiclass_stranded_STAR', args.psiclass_stranded_STAR_path),
            ('stringtie_stranded_default_STAR', args.stringtie_stranded_default_STAR_path),
            ('stringtie_stranded_AltCommands_STAR', args.stringtie_stranded_AltCommands_STAR_path)
        ]
    }

    # Add optional evidences
    optional_evidences = {
        "stringtie_Isoseq_default": args.stringtie_Isoseq_default_path,
        "stringtie_Isoseq_AltCommands": args.stringtie_Isoseq_AltCommands_path,
        "psiclass_unstranded_STAR": args.psiclass_unstranded_STAR_path,
        "stringtie_unstranded_default_STAR": args.stringtie_unstranded_default_STAR_path,
        "stringtie_unstranded_AltCommands_STAR": args.stringtie_unstranded_AltCommands_STAR_path,
    }

    for key, path in optional_evidences.items():
        if path:
            evidences['transcriptome'].append((key, path))

    # Print full evidences dictionary
    print("-------------------- Complete evidences used:")
    for category, items in evidences.items():
        print(f"  {category}:")
        for name, path in items:
            print(f"    {name}: {path}")

    if chosen_chromosomes:
        # Create annotation objects
        print("-------------------- Creating ab initio annotations")
        ab_initio_annotations = [Annotation(name, path, genome, chosen_chromosomes=chosen_chromosomes) for name, path in evidences['ab_initio']]

        # Create annotation objects
        print("-------------------- Creating transcriptome annotations")
        transcriptome_annotations = [Annotation(name, path, genome, chosen_chromosomes=chosen_chromosomes) for name, path in evidences['transcriptome']]
    else:
        # Create annotation objects
        print("-------------------- Creating ab initio annotations")
        ab_initio_annotations = [Annotation(name, path, genome) for name, path in evidences['ab_initio']]

        # Create annotation objects
        print("-------------------- Creating transcriptome annotations")
        transcriptome_annotations = [Annotation(name, path, genome) for name, path in evidences['transcriptome']]

    # Merge evidences
    print("-------------------- Merging annotations")
    merged_annotation = transcriptome_annotations[0].copy()
    for annotation in transcriptome_annotations[1:]:
        print(f"Merging transcriptome annotation: {annotation.name}")
        merged_annotation.merge(annotation)

    for annotation in ab_initio_annotations:
        print(f"Merging ab initio annotation: {annotation.name}")
        merged_annotation.merge(annotation)

    print("-------------------- Processing alternative transcripts into genes")
    merged_annotation.make_alternative_transcripts_into_genes()

    if chosen_chromosomes:
        merged_annotation.id = f"merged_annotation_{'_'.join(chosen_chromosomes)}"
        merged_annotation.name = f"merged_annotation_transcripts_{'_'.join(chosen_chromosomes)}"
    else:
        merged_annotation.id = "merged_annotation"
        merged_annotation.name = 'merge_annotation_transcripts'

    print(f"-------------------- Exporting merged GFF to {args.output_dir}")
    # merged_annotation.export_gff(args.output_dir, args.output_gff)
    merged_annotation.export_gff(args.output_dir)

    print("-------------------- Exporting merged unique proteins")
    merged_annotation.export_unique_proteins(custom_path=args.output_dir, genome=genome)

    print("-------------------- Updating stats and exporting")
    merged_annotation.update_stats(export=True, genome=genome)

    # Save pickle (to load the merged object with all annotations in case anything fails in the following steps)
    print(f"-------------------- Saving merged annotation to {args.output_pickle}")
    pickle_save(args.output_pickle, merged_annotation)

# Parse arguments
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Genome annotation merging.")
    parser.add_argument("--genome_name", required=True, help="Name of the genome - MANDATORY")
    parser.add_argument("--genome_path", required=True, help="Path to the genome FASTA file - MANDATORY")
    parser.add_argument("--chromosome", required=False, help="Chromosome ID to run Aegis per chromosome - OPTIONAL - if not given, all the chromosome are treated together.")
    parser.add_argument("--augustus_path", required=True, help="Path to Augustus GFF3 file - MANDATORY")
    parser.add_argument("--genemark_path", required=True, help="Path to GeneMark GFF3 file - MANDATORY")
    parser.add_argument("--liftoff_path", required=True, help="Path to Liftoff GFF3 file - MANDATORY")
    parser.add_argument("--psiclass_stranded_STAR_path", required=True, help="Path to Psiclass stranded GFF3 file - MANDATORY")
    parser.add_argument("--stringtie_stranded_default_STAR_path", required=True, help="Path to StringTie stranded default GFF3 file - MANDATORY")
    parser.add_argument("--stringtie_stranded_AltCommands_STAR_path", required=True, help="Path to StringTie stranded AltCommands GFF3 file - MANDATORY")
    parser.add_argument("--psiclass_unstranded_STAR_path", required=False, help="Path to Psiclass unstranded GFF3 file - OPTIONAL")
    parser.add_argument("--stringtie_unstranded_default_STAR_path", required=False, help="Path to StringTie unstranded default GFF3 file - OPTIONAL")
    parser.add_argument("--stringtie_unstranded_AltCommands_STAR_path", required=False, help="Path to StringTie unstranded AltCommands GFF3 file - OPTIONAL")
    parser.add_argument("--stringtie_Isoseq_default_path", required=False, help="Path to StringTie Isoseq default GFF3 file - OPTIONAL")
    parser.add_argument("--stringtie_Isoseq_AltCommands_path", required=False, help="Path to StringTie Isoseq AltCommands GFF3 file - OPTIONAL")
    parser.add_argument("--output_dir", required=True, help="Directory to save the outputs - MANDATORY")
    parser.add_argument("--output_gff", required=True, help="Name of the merged GFF file - MANDATORY")
    parser.add_argument("--output_pickle", required=True, help="Path to save the merged pickle file - MANDATORY")
    args = parser.parse_args()

    main(args)

# example of complete bash command to launch the script
# python3 script.py \
#   --genome_name PN40024_T2T \
#   --genome_path /storage/tom/PN40024_T2T.fasta \
#   --chromosome chr01 \
#   --augustus_path /storage/tom/augustus_T2T/Augustus/augustus.hints.gff3 \
#   --genemark_path /storage/tom/augustus_T2T/GeneMark-ETP/genemark.gff3 \
#   --liftoff_path /storage/tom/inigo/annotation/T2T/liftoff/soft_masked/T2T_ref_soft_masked.gff3 \
#   --psiclass_stranded_STAR_path /storage/tom/T2T_mapping_bam/merged_T2T_PN40024_stranded_psiclass.gff3 \
#   --stringtie_stranded_default_STAR_path /storage/tom/T2T_mapping_bam/merged_T2T_PN40024_stranded_stringtie_default.gff3 \
#   --stringtie_stranded_AltCommands_STAR_path /storage/tom/T2T_mapping_bam/merged_T2T_PN40024_stranded_stringtie_morus.gff3 \
#   --psiclass_unstranded_STAR_path /storage/tom/T2T_mapping_bam/merged_T2T_PN40024_unstranded_psiclass.gff3 \
#   --stringtie_unstranded_default_STAR_path /storage/tom/T2T_mapping_bam/merged_T2T_PN40024_unstranded_stringtie_default.gff3 \
#   --stringtie_unstranded_AltCommands_STAR_path /storage/tom/T2T_mapping_bam/merged_T2T_PN40024_unstranded_stringtie_morus.gff3 \
#   --stringtie_Isoseq_default_path /storage/tom/T2T_mapping_bam/Isoseq/hq_transcripts.RI_rmv_on_T2T_sorted_stringtieDefaultCommands.gff3 \
#   --stringtie_Isoseq_AltCommands_path /storage/tom/T2T_mapping_bam/Isoseq/hq_transcripts.RI_rmv_on_T2T_sorted_stringtieAltCommands.gff3 \
#   --output_dir /storage/tom/chr01 \
#   --output_gff merged_annotation_chr01.gff3 \
#   --output_pickle /storage/tom/merged_annotation_chr01.pkl