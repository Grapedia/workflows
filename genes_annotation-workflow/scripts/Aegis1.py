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
    # Upload the genome assembly
    print(f"Loading genome: {args.genome_name}")
    genome = Genome(args.genome_name, args.genome_path)

    # define evidences
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

    # add stringtie_Isoseq_default_path if given
    if args.stringtie_Isoseq_default_path:
        evidences['transcriptome'].append(('stringtie_Isoseq_default', args.stringtie_Isoseq_default_path))

    # add stringtie_Isoseq_AltCommands_path if given
    if args.stringtie_Isoseq_AltCommands_path:
        evidences['transcriptome'].append(('stringtie_Isoseq_AltCommands', args.stringtie_Isoseq_AltCommands_path))

    # add psiclass_unstranded_STAR_path if given
    if args.psiclass_unstranded_STAR_path:
        evidences['transcriptome'].append(('psiclass_unstranded_STAR', args.psiclass_unstranded_STAR_path))

    # add stringtie_unstranded_default_STAR_path if given
    if args.stringtie_unstranded_default_STAR_path:
        evidences['transcriptome'].append(('stringtie_unstranded_default_STAR', args.stringtie_unstranded_default_STAR_path))

    # add stringtie_unstranded_AltCommands_STAR_path if given
    if args.stringtie_unstranded_AltCommands_STAR_path:
        evidences['transcriptome'].append(('stringtie_unstranded_AltCommands_STAR', args.stringtie_unstranded_AltCommands_STAR_path))

    # Create annotation objects
    print("Creating ab initio annotations")
    ab_initio_annotations = [Annotation(name, path, genome) for name, path in evidences['ab_initio']]

    # Create annotation objects
    print("Creating transcriptome annotations")
    transcriptome_annotations = [Annotation(name, path, genome) for name, path in evidences['transcriptome']]

    # Merge evidences
    print("Merging annotations")
    merged_annotation = transcriptome_annotations[0].copy()
    for annotation in transcriptome_annotations[1:]:
        merged_annotation.merge(annotation)

    for annotation in ab_initio_annotations:
        merged_annotation.merge(annotation)

    if args.annotation_pickle:
        print(f"Loading pickle file: {args.annotation_pickle}")
        merged_annotation = pickle_load(args.annotation_pickle)

    print("Processing alternative transcripts into genes")
    merged_annotation.make_alternative_transcripts_into_genes()

    merged_annotation.id = 'merge_annotation'
    merged_annotation.name = 'merge_annotation_transcripts'

    print(f"Exporting merged GFF to {args.output_dir}")
    merged_annotation.export_gff(args.output_dir, args.output_gff)

    print("Exporting merged unique proteins")
    merged_annotation.export_unique_proteins(custom_path=args.output_dir, genome=genome)

    print("Updating stats and exporting")
    merged_annotation.update_stats(export=True, genome=genome)

    # Save pickle (to load the merged object with all annotations in case anything fails in the following steps)
    print(f"Saving merged annotation to {args.output_pickle}")
    pickle_save(args.output_pickle, merged_annotation)

# Parse arguments
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Genome annotation merging.")
    parser.add_argument("--genome_name", required=True, help="Name of the genome - MANDATORY")
    parser.add_argument("--genome_path", required=True, help="Path to the genome FASTA file - MANDATORY")
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
    parser.add_argument("--annotation_pickle", required=False, help="Path to an existing annotation pickle file")
    parser.add_argument("--output_dir", required=True, help="Directory to save the outputs")
    parser.add_argument("--output_gff", required=True, help="Name of the merged GFF file")
    parser.add_argument("--output_pickle", required=True, help="Path to save the merged pickle file")
    args = parser.parse_args()

    main(args)

# example of complete bash command to launch the script
# python3 script.py \
#   --genome_name PN40024_T2T \
#   --genome_path /storage/tom/PN40024_T2T.fasta \
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
#   --output_dir /storage/tom/ \
#   --output_gff merged_annotation.gff3 \
#   --output_pickle /storage/tom/merged_annotation.pkl