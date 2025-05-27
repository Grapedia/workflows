#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import time
import os
from modules.tools import pickle_load, pickle_save
from modules.annotation import Annotation
from modules.genome import Genome

def main(args):

    chosen_chromosomes = [args.chromosome] if args.chromosome else None
    if chosen_chromosomes:
        print(f"-------------------- Chosen chromosome is : {chosen_chromosomes}")
    else:
        print("-------------------- Aegis launched on all chromosomes.")

    # Load the merged annotation pickle
    print('Loading merged annotation pickle file ...')
    merge = pickle_load(args.merged_annotation)

    print('Load hard-masked genome...')
    assembly_hard_masked = Genome('assembly_hard_masked', args.hard_masked_genome)
    merge.calculate_transcript_masking(hard_masked_genome=assembly_hard_masked)

    if args.update:
        print('## Probably not needed in the future ## - Update object to sort it')
        merge.update()

    print('Overlaps...')
    # merge.detect_gene_overlaps()
    merge.detect_gene_overlaps(sort_processes=1)


    print('Adding DIAMOND results...')
    diamond_hits_dict = {}
    for item in args.diamond_hits:
        name, path = item.split("=")
        diamond_hits_dict[name] = path

    for name, path in diamond_hits_dict.items():
        print(f"Adding of {name} results from {path}")
        merge.add_blast_hits(name, path)

    # Save pickle
    print('Saving intermediate annotation...')
    pickle_save(args.intermediate_annotation, merge)

    print('Exporting GFF...')
    merge.update()

    if chosen_chromosomes:
        merged_annotation_blast_filename = f"merged_annotation_blast_{'_'.join(chosen_chromosomes)}"
        merge.export_gff(args.export_dir, merged_annotation_blast_filename)
    else:
        merge.export_gff(args.export_dir, 'merged_annotation_blast.gff3')

    print('Reducing redundancy...')
    merge.remove_redundancy(source_priority=args.source_priority, hard_masked_genome=assembly_hard_masked)

    # Save final pickle
    print('Saving final annotation...')
    pickle_save(args.final_annotation, merge)
    if chosen_chromosomes:
        merge.id = f"final_annotation_{'_'.join(chosen_chromosomes)}"
        merge.name = f"final_annotation_{'_'.join(chosen_chromosomes)}"
    else:
        merge.id = 'final_annotation'
        merge.name = 'final_annotation'

    print('Final exports...')
    if chosen_chromosomes:
        final_annotation_filename = f"final_annotation_{'_'.join(chosen_chromosomes)}.gff3"
        merge.export_gff(args.final_export_dir, final_annotation_filename)
        merge.export_equivalences(stringent=False, return_df=False, export_csv=True, export_self=True, NAs=False)
        merge.generate_sequences(genome=assembly_hard_masked, just_CDSs=True)
        merge.export_proteins(only_main=True, verbose=False, custom_path=args.final_export_dir)
        merge.export_proteins(only_main=False, verbose=False, custom_path=args.final_export_dir)
    else:
        merge.export_gff(args.final_export_dir, 'final_annotation.gff3')
        merge.export_equivalences(stringent=False, return_df=False, export_csv=True, export_self=True, NAs=False)
        merge.generate_sequences(genome=assembly_hard_masked, just_CDSs=True)
        merge.export_proteins(only_main=True, verbose=False, custom_path=args.final_export_dir)
        merge.export_proteins(only_main=False, verbose=False, custom_path=args.final_export_dir)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Automate genome annotation pipeline.")

    parser.add_argument('--merged_annotation', type=str, required=True, help="Path to the merged_annotation.pkl file.")
    parser.add_argument('--hard_masked_genome', type=str, required=True, help="Path to the hard-masked genome file.")
    parser.add_argument('--diamond_hits', nargs='+', required=True, help="List of DIAMOND results in the format 'name=path'. Example : --diamond_hits Viridiplantae=/path/to/viridiplantae_vs_proteins_assembly.diamond Eudicots=/path/to/eudicots_vs_proteins_assembly.diamond Araport=/path/to/araport_vs_proteins_assembly.diamond")
    parser.add_argument('--intermediate_annotation', type=str, required=True, help="Path to save the intermediate annotation pickle file.")
    parser.add_argument('--final_annotation', type=str, required=True, help="Path to save the final annotation pickle file.")
    parser.add_argument('--export_dir', type=str, required=True, help="Directory to export intermediate GFF files.")
    parser.add_argument('--final_export_dir', type=str, required=True, help="Directory to export final GFF and other outputs.")
    parser.add_argument('--update', action='store_true', help="Update the object during the pipeline.")
    parser.add_argument("--chromosome", required=False, help="Chromosome ID to run Aegis per chromosome - OPTIONAL - if not given, all the chromosome are treated together.")
    parser.add_argument('--source_priority', nargs='+', default=['Araport', 'Viridiplantae', 'Eudicots'], help="Source priority list for redundancy reduction.")

    args = parser.parse_args()

    main(args)
