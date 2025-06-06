"""
Created on Fri Oct 7 15:09:46 2022

Module defining several genomic classes.

@authors: David Navarro, Antonio Santiago
"""

import sys
import copy
import random
import time
import pandas as pd
import matplotlib.pyplot as plt
#import networkx as nx

from modules.plots import pie_chart, barplot
from modules.genefunctions import overlap, reverse_complement, find_all_occurrences
from modules.feature import Feature
from modules.gene import Gene
from modules.transcript import Transcript
from modules.subfeatures import Exon, UTR, Intron


from modules.hits import OverlapHit, BlastHit

from os import system
from statistics import mean
from scipy.stats import fisher_exact
from tqdm import tqdm
from multiprocessing import Pool

def sort_and_update_genes(chrom, genes_dict):
    genes = sorted(genes_dict.values())
    sorted_genes = {g.id: g.copy() for g in genes}
    return chrom, sorted_genes

class Annotation():
    # These categories must be updated based on different gff terms, since 
    # these are not fixed, the more gffs are looked at the more terms appear
    # so add them here to be read into the annotation object as long as they
    # are unambiguous
    features = {}
    features["gene"] = ["gene", "pseudogene", "transposable_element_gene"]
    coding_transcripts = ["mRNA"]
    noncoding_transcripts = ["antisense_lncRNA", "antisense_RNA", 
                             "miRNA_primary_transcript", "ncRNA", "lncRNA",
                             "lnc_RNA", "pseudogenic_tRNA", "rRNA", "snoRNA",
                             "snRNA", "tRNA"]
    # Some features are clearly transcript level features but they cannot be
    # classed as coding/noncoding just by looking at the name   
    features["transcript"] = (["transcript", "transcript_region", "primary_transcript",
                               "pseudotranscript", "pseudogenic_transcript", "mRNA_TE_gene"] 
                               + coding_transcripts + noncoding_transcripts)
    features["UTR"] = ["UTR", "three_prime_UTR", "five_prime_UTR", "five_prime_utr", "three_prime_utr"]
    features["exon"] = ["exon", "pseudogenic_exon"]
    features["CDS"] = ["CDS"]
    features["other_subfeature"] = ["miRNA"]
    features_r = {}
    for key, values in features.items():
        for value in values:
            features_r[value] = key
    bar_colors = ["31", "32", "33", "33", "33", "33", "34"]
    def __init__(self, name:str, annot_file_path:str, genome:object, original_annotation:object=None, target:bool=False, to_overlap:bool=True, rework_CDSs:bool=True, chosen_chromosomes:list=None, chosen_coordinates:tuple=None, sort_processes:int=2, define_synteny=False):
        
        start = time.time()
        if chosen_chromosomes != None:
            if len(chosen_chromosomes) > 1:
                self.name = f"{name}_{chosen_chromosomes[0]}-{chosen_chromosomes[-1]}"
            else:
                self.name = f"{name}_{chosen_chromosomes[0]}"
        else:
            self.name = name
        self.file = annot_file_path
        self.path = "/".join(self.file.split("/")[0:-1]) + "/"
        self.gff_header = ""
        self.target = target
        self.to_overlap = to_overlap
        self.liftoff = False
        self.merged = False
        self.sorted = False
        self.genome = genome.name
        self.id = f"{self.name}_on_{self.genome}"

        print(f"\nProcessing {self.id} annotation object\n")
        self.excluded_chromosomes = []
        self.features = set()
        # genes will be added as {"ch":{"gene_id" : gene_object}}
        self.chrs = {}
        # we save here chr as the value corresponding to each gene id
        self.all_gene_ids = {}
        # we save here chr as the value corresponding to each transcript id
        self.all_transcript_ids = {}
        self.all_protein_ids = {}
        self.unmapped = []
        # Here we insert any feature which did not fit the current standards
        self.atypical_features = []
        self.self_overlapping = []
        self.overlapped_annotations = []

        self.feature_suffix = ""
        self.suffix = ""

        self.coding_removed = False
        if "_minus_coding" in annot_file_path:
            self.coding_removed = True

        self.non_coding_removed = False
        if "_minus_non_coding" in annot_file_path:
            self.non_coding_removed = True

        self.transposable_removed = False
        if "_minus_TE" in annot_file_path:
            self.transposable_removed = True

        self.non_transposable_removed = False
        if "_minus_non_TE" in annot_file_path:
            self.non_transposable_removed = True

        self.featurecounts = False
        if "_fcounts" in annot_file_path:
            self.featurecounts = True

        self.aegis = False
        if "_aegis" in annot_file_path:
            self.aegis = True

        self.combined = False
        if "_combined" in annot_file_path:
            self.combined = True

        self.clean = False
        if "_clean" in annot_file_path:
            self.clean = True

        self.dapmod = False
        if "_dapmod" in annot_file_path:
            self.dapmod = True

        self.confrenamed = False
        if "_confrenamed" in annot_file_path:
            self.confrenamed = True

        self.dapfit = genome.dapfit

        self.symbols_added = False

        self.renamed_ids = False

        self.small_cds_removed = False
        if "_minus_small_CDSs" in annot_file_path:
            self.small_cds_removed = True

        self.promoter_size = 3000
        
        self.generated_all_sequences = False
        self.contains_all_sequences = False        
        self.generated_CDS_sequences = False
        self.contains_CDS_sequences = False
        self.generated_protein_sequences = False
        self.contains_protein_sequences = False

        self.promoter_types = "standard"
        if original_annotation != None:
            self.liftoff = True

        self.errors = {"repeated_gene_IDs": [], "1bp_gene": [],
                       "1bp_transcript": [], "subfeature_to_gene": [],
                       "transcript_to_more_than_1_gene": [],
                       "transcript_to_inexistent_gene": [],
                       "repeat_transcript_different_genes": [],
                       "repeat_transcript_same_gene": [],
                       "transcript_to_gene_other_chr": []}
        
        self.warnings = {"1bp_exon": [], "1bp_CDS": [], "1bp_UTR": [], "decreasing_coordinates": [],
                         "missing_subfeature_parent": [],
                         "missing_subfeature_parent_liftoff": [],
                         "multiple_CDSs_per_transcript": [],
                         "possible_policistronic_transcript": [],
                         "transcript_with_no_exons": [],
                         "gene_with_no_transcripts": []}
        
        parsed_lines = {key: [] for key in Annotation.features}
        parsed_lines["atypical"] = []
        
        temp = []
        protein_match_lines = []
        chromosomes_t = set()
        # Read lines from input file, filtering out lines consisting of "###\n"
        with open(self.file, encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if line == "":
                    continue
                if line.startswith("#") and line != "###" and line != "#":
                    self.gff_header += line + "\n"
                elif line != "###" and line != "#":
                    line = line.split("\t")
                    if len(line) > 2:
                        if chosen_chromosomes != None:
                            if line[0] not in chosen_chromosomes:
                                continue
                            else:
                                if chosen_coordinates != None:
                                    if int(line[3]) < chosen_coordinates[0] or int(line[4]) > chosen_coordinates[1]:
                                        continue
                        self.features.add(line[2])
                        if line[2] == "nucleotide_to_protein_match":
                            #to check why ":"
                            chromosomes_t.add(line[0].split(":")[0])
                            protein_match_lines.append(line)
                        elif line[2] in Annotation.features_r:
                            chromosomes_t.add(line[0])
                            temp.append(line)
                        else:
                            chromosomes_t.add(line[0])
                            parsed_lines["atypical"].append(line)
                    else:
                        self.gff_header += "#" + "\t".join(line) + "\n"

        sorted_lines = sorted(temp, key=lambda x: (int(x[3]), int(x[4])))
        protein_match_lines = sorted(protein_match_lines, key=lambda x: (int(x[3]), int(x[4])))
        temp.clear()
        for line in sorted_lines:
            parsed_lines[Annotation.features_r[line[2]]].append(line)
        sorted_lines.clear()
        chromosomes_t = list(chromosomes_t)
        chromosomes_t.sort()
        for ch in chromosomes_t:
            self.chrs[ch] = {}
        self.features = list(self.features)

        

        # Check if stdout or stderr are redirected to files
        stdout_redirected = not sys.stdout.isatty()
        stderr_redirected = not sys.stderr.isatty()

        # Disable tqdm if stdout or stderr are redirected
        if stdout_redirected or stderr_redirected:
            disable = True
        else:
            disable = False

        created_genes = False
        transcript_info = {}

        if self.features == ["nucleotide_to_protein_match"]:
            print(f"{self.id} annotation has just nucleotide to protein matches so genes and transcripts will be generated")
            progress_bar = tqdm(total=len(protein_match_lines),
                    disable=disable,
                    bar_format=(
            f'\033[1;{Annotation.bar_colors[0]}mReading protein match features:\033[0m '
            '{percentage:3.0f}%|'
            f'\033[1;{Annotation.bar_colors[0]}m{{bar}}\033[0m| '
            '{n}/{total} [{elapsed}<{remaining}]'))
            temp_protein_matches = {}
            for line in protein_match_lines:
                progress_bar.update(1)
                ch = line[0].split(":")[0]
                source = line[1]
                ft = line[2]
                coord = [int(line[3]), int(line[4])]
                if coord[0] > coord[1]:
                    print(f"{self.id} Warning: Decreasing coordinates for {ft} {ID}")
                    coord.sort()
                score = line[5]
                strand = line[6]
                phase = line[7]
                attributes = line[8].strip(";")
                attributes_l = attributes.split(";")
                ID = False
                for a in attributes_l:
                    a = a.split("=")
                    if a[0] == "ID":
                        ID = a[1].strip()
                if not ID:
                    print(f"Error: ID not found for {ft}: {line}")
                    continue
                if ID not in temp_protein_matches:
                    temp_protein_matches[ID] = [Feature(ID, ch, source, ft, strand,coord[0], coord[1], score, phase, attributes)]
                else:
                    temp_protein_matches[ID].append(Feature(ID, ch, source, ft, strand,coord[0], coord[1], score, phase, attributes))
            
            progress_bar.close()
            progress_bar = tqdm(total=len(temp_protein_matches.keys()),
                    disable=disable,
                    bar_format=(
            f'\033[1;{Annotation.bar_colors[0]}mCreating genes for nucleotide to protein matches:\033[0m '
            '{percentage:3.0f}%|'
            f'\033[1;{Annotation.bar_colors[0]}m{{bar}}\033[0m| '
            '{n}/{total} [{elapsed}<{remaining}]'))

            for ID, features in temp_protein_matches.items():
                progress_bar.update(1)
                count = 0
                gene_id = f"{ID}_gene"
                sorted_features = features.copy()
                sorted_features.sort()
                ch_size = genome.scaffolds[sorted_features[0].ch].size
                if gene_id not in self.all_gene_ids:
                    self.chrs[sorted_features[0].ch][gene_id] = Gene(False, False,
                                                        gene_id, sorted_features[0].ch, ch_size,
                                                        sorted_features[0].source, "gene",
                                                        sorted_features[0].strand, sorted_features[0].start, sorted_features[-1].end,
                                                        sorted_features[0].score, ".", f"ID={gene_id}")
                    self.all_gene_ids[gene_id] = sorted_features[0].ch
                else:
                    while gene_id in self.all_gene_ids:
                        count += 1
                        gene_id = f"{ID}_{count}_gene"
                    self.chrs[sorted_features[0].ch][gene_id] = Gene(False, False,
                                                        gene_id, sorted_features[0].ch, ch_size,
                                                        sorted_features[0].source, "gene",
                                                        sorted_features[0].strand, sorted_features[0].start, sorted_features[-1].end,
                                                        sorted_features[0].score, ".", f"ID={gene_id}")
                    self.all_gene_ids[gene_id] = sorted_features[0].ch
                
                self.chrs[sorted_features[0].ch][gene_id].transcripts["temp_t"] = Transcript("temp_t", sorted_features[0].ch, ch_size, sorted_features[0].source, "mRNA",
                                                        sorted_features[0].strand, sorted_features[0].start, sorted_features[-1].end,
                                                        sorted_features[0].score, ".", f"ID=temp_t;Parent={gene_id}")
                self.all_transcript_ids["temp_t"] = sorted_features[0].ch
                for f in sorted_features:
                    self.chrs[f.ch][gene_id].transcripts["temp_t"].temp_CDSs.append(Feature("temp_CDS", f.ch, f.source, "CDS", f.strand, f.start, f.end, f.score, f.phase, f"ID=temp_CDS;Parent=temp_t"))
                    self.chrs[f.ch][gene_id].transcripts["temp_t"].exons.append(Feature(f"temp_exon", f.ch, f.source, "exon", f.strand, f.start, f.end, f.score,  ".", f"ID=temp_exon;Parent=temp_t"))

            progress_bar.close()

        else:
            if "gene" not in self.features:
                created_genes = True
                progress_bar = tqdm(total=len(parsed_lines["transcript"]),
                                    disable=disable,
                                    bar_format=(
                        f'\033[1;{Annotation.bar_colors[0]}mCreating Genes:\033[0m '
                        '{percentage:3.0f}%|'
                        f'\033[1;{Annotation.bar_colors[0]}m{{bar}}\033[0m| '
                        '{n}/{total} [{elapsed}<{remaining}]'))
                for line in parsed_lines["transcript"]:
                    progress_bar.update(1)
                    ch = line[0]
                    source = line[1]
                    ft = line[2]
                    coord = [int(line[3]), int(line[4])]
                    score = line[5]
                    strand = line[6]
                    phase = line[7]
                    attributes = line[8].strip(";")
                    attributes_l = attributes.split(";")
                    ch_size = genome.scaffolds[ch].size
                    ID = False
                    for a in attributes_l:
                        a = a.split("=")
                        if a[0] == "ID":
                            ID = a[1].strip()
                    if not ID:
                        print(f"Error: ID not found for {ft}: {line}")
                        continue

                    if "pseudo" in ft:
                        pseudogene = True
                    else:
                        pseudogene = False

                    if ft == "transposable_element_gene":
                            transposable = True
                    elif "transposable=True" in attributes:
                        transposable = True
                    else:
                        transposable = False                  

                    #improve at some point, slow with massive lists
                    if ID not in self.all_transcript_ids:
                        self.all_transcript_ids[f"{ID}"] = ch
                        self.all_gene_ids[f"{ID}_gene"] = ch
                        transcript_info[f"{ID}"] = [ch, coord[0], coord[1], 1]
                        self.chrs[ch][f"{ID}_gene"] = Gene(pseudogene, transposable,
                                                    f"{ID}_gene", ch, ch_size,
                                                    source, "gene",
                                                    strand, coord[0], coord[1],
                                                    score, ".", f"ID={ID}_gene")
                    else:
                        transcript_info[f"{ID}"][3] += 1
                        self.all_gene_ids[f"{ID}_{transcript_info[ID][3]}_gene"] = ch
                        self.chrs[ch][f"{ID}_{transcript_info[ID][3]}_gene"] = Gene(pseudogene, transposable,
                                                f"{ID}_{transcript_info[ID][3]}_gene",
                                                ch, ch_size, source, "gene", 
                                                strand, coord[0], coord[1], score,
                                                 ".",
                                                f"ID={ID}_{transcript_info[ID][3]}_gene")
                        self.all_transcript_ids[f"{ID}_{transcript_info[ID][3]}"] = ch
                        print(f"{self.id} Warning: repeated transcript ID "
                            f"{ID}_gene renamed to {ID}_{transcript_info[ID][3]}_gene")

                progress_bar.close()

                print(f"{self.id} did not have genes and these were created from transcripts")

            for n, (ft_level, lines) in enumerate(parsed_lines.items()):
                progress_bar = tqdm(total=len(lines), disable=disable, bar_format=(
                            f'\033[1;{Annotation.bar_colors[n]}mAdding {ft_level}s:\033[0m '
                            '{percentage:3.0f}%|'
                            f'\033[1;{Annotation.bar_colors[n]}m{{bar}}\033[0m| '
                            '{n}/{total} [{elapsed}<{remaining}]'))
                for line in lines:
                    ch = line[0]
                    source = line[1]
                    ft = line[2]
                    coord = [int(line[3]), int(line[4])]
                    score = line[5]
                    strand = line[6]
                    phase = line[7]
                    attributes = line[8].strip(";")
                    attributes_l = attributes.split(";")
                    ID = False
                    for a in attributes_l:
                        a = a.split("=")
                        if a[0] == "ID":
                            ID = a[1].strip()
                    if not ID:
                        ID = ""

                    if coord[0] > coord[1]:
                        print(f"{self.id} Warning: Decreasing coordinates "
                                f"for {ft} {ID}")
                        self.warnings["decreasing_coordinates"].append(ID)
                        coord.sort()

                    # genes
                    if n == 0:
                        ch_size = genome.scaffolds[ch].size
                        if coord[0] == coord[1]:
                            print(f"{self.id} Error: 1bp gene {ID} feature")
                            self.errors["1bp_gene"].append(ID)
                        if ft == "pseudogene":
                            pseudogene = True
                        else:
                            pseudogene = False

                        if ft == "transposable_element_gene":
                            transposable = True
                        elif "transposable=True" in attributes:
                            transposable = True
                        else:
                            transposable = False

                        #improve at some point, slow with massive lists
                        if ID not in self.all_gene_ids:
                            self.all_gene_ids[ID] = ch
                            self.chrs[ch][ID] = Gene(pseudogene, transposable,
                                                    ID, ch, ch_size,
                                                    source, ft, strand, coord[0],
                                                    coord[1], score, ".",
                                                    attributes)
                        else:
                            print(f"{self.id} Error: repeated gene ID {ID}")
                            self.errors["repeated_gene_IDs"].append(ID)

                    # transcripts with created genes
                    elif n == 1 and created_genes:
                        ch_size = genome.scaffolds[ch].size
                        if coord[0] == coord[1]:
                            print(f"{self.id} Error: 1bp transcript level {ID} "
                                "feature ")
                            self.errors["1bp_transcript"].append(ID)  
                        count =  transcript_info[f"{ID}"][3]

                        if count == 1:
                            self.chrs[ch][f"{ID}_gene"].transcripts[ID] = Transcript(ID, ch, ch_size, source, ft, 
                                                            strand, coord[0],
                                                            coord[1], score,
                                                            ".",
                                                        f"ID={ID};Parent={ID}_gene")
                        else:
                            found = False
                            for x in range(count):
                                if x == 0:
                                    t_parent = f"{ID}_gene"
                                    t_id = f"{ID}"
                                else:
                                    t_parent = f"{ID}_{x+1}_gene"
                                    t_id = f"{ID}_{x+1}"

                                if t_parent in self.chrs[ch]:
                                    start_parent = self.chrs[ch][t_parent].start
                                    end_parent = self.chrs[ch][t_parent].end
                                    if ((start_parent <= coord[0])
                                        and (end_parent >= coord[1])):
                                        found = True
                                        (self.chrs[ch][t_parent]
                                        .transcripts[t_id]) = Transcript(t_id, ch,
                                                                        ch_size,
                                                                    source, ft,
                                                                    strand,
                                                                    coord[0],
                                                                    coord[1], score,
                                                                    ".",
                                                        f"ID={t_id};Parent={t_parent}")
                                        break
                            if not found:
                                print("Error: repeated transcript ID {ID} could "
                                    "not find its newly created gene copy")
    
                    # transcripts
                    elif n == 1:
                        ch_size = genome.scaffolds[ch].size
                        if coord[0] == coord[1]:
                            print(f"{self.id} Error: 1bp transcript level {ID} "
                                "feature ")
                            self.errors["1bp_transcript"].append(ID)
                        for term in attributes_l:
                            if "arent=" not in term:
                                continue # skip term if not parent term
                            term = term.split("=")[1]
                            parents = term.split(",")
                            if len(parents) > 1: # does not happen it seems
                                print(f"{self.id} Error: {ID} transcript refers "
                                        "to more than 1 gene!")
                                (self.errors["transcript_to_more_than_1_gene"]
                                        .append(ID))
                            else:
                                parent = parents[0].strip()
                                if parent in self.chrs[ch]:
                                    if ID in self.chrs[ch][parent].transcripts:
                                        print(f"{self.id} Error: duplicated {ID}"
                                                f"transcript for {parent} gene")
                                        (self.errors["repeat_transcript_same_gene"]
                                                    .append(ID))
                                    else:
                                        (self.chrs[ch][parent].transcripts
                                                [ID]) = Transcript(ID, ch, ch_size, source,
                                                                ft, 
                                                                strand, coord[0],
                                                                coord[1], score,
                                                                ".",
                                                                attributes)
                                        self.all_transcript_ids[ID] = ch
                                elif parent in self.all_gene_ids:
                                    print(f"{self.id} Error: {ID} transcript " 
                                            "refers to a gene in a different "
                                            "chromosome")
                                    (self.errors["transcript_to_gene_other_chr"]
                                            .append(ID))
                                else:
                                    print(f"{self.id} Error: {ID} transcript " 
                                            "refers to an inexistent gene")
                                    (self.errors["transcript_to_inexistent_gene"]
                                            .append(ID))
                                    
                    # transcript subfeatures
                    elif n < (len(parsed_lines) - 1):
                        if coord[0] == coord[1] and n == 3:
                            self.warnings["1bp_exon"].append(ID)
                        if coord[0] == coord[1] and n == 4:
                            self.warnings["1bp_CDS"].append(ID)
                        if coord[0] == coord[1] and n == 2:
                            self.warnings["1bp_UTR"].append(ID)

                        for term in attributes_l:
                            if "arent=" not in term:
                                continue # skip term if not parent term
                            term = term.split("=")[1]
                            parents = term.split(",")
                            for parent in parents:
                                if not parent:
                                    continue
                                parent = parent.strip()
                                #improve at some point, slow with massive lists
                                # parent not found within transcripts
                                if parent not in self.all_transcript_ids:
                                    if parent in self.chrs[ch]:
                                        # subset of cases where subfeature
                                        # directly references a pseudogene
                                        if self.chrs[ch][parent].pseudogene:
                                            # creating a pseudotranscript
                                            # for the pseudogene
                                            if parent[0:4] == "gene":
                                                pseudo_t = f"pseudo_t_{parent[5:]}"
                                            else:
                                                pseudo_t = f"pseudo_t_{parent}"

                                            if self.chrs[ch][parent].transcripts == {}:
                                                # PSEUDOGENES refered to by
                                                # transcript subfeatures
                                                # are given a single
                                                # pseudotranscript
                                                self.chrs[ch][parent].transcripts[pseudo_t] = Transcript(pseudo_t, self.chrs[ch][parent].ch, self.chrs[ch][parent].ch_size, self.chrs[ch][parent].source, "pseudotranscript", self.chrs[ch][parent].strand, self.chrs[ch][parent].start, self.chrs[ch][parent].end, self.chrs[ch][parent].score, ".", self.chrs[ch][parent].attributes)
                                                
                                            if pseudo_t in self.chrs[ch][parent].transcripts:
                                                if n == 4:
                                                    # CDS segments
                                                    (self.chrs[ch][parent]
                                                    .transcripts[pseudo_t]
                                                    .temp_CDSs.append(Feature
                                                            (ID, ch, source, ft, 
                                                            strand, coord[0],
                                                            coord[1], score, phase, 
                                                            attributes)))
                                                elif n == 3:
                                                    (self.chrs[ch][parent]
                                                    .transcripts[pseudo_t]
                                                    .exons.append
                                                    (Exon(ID, ch, source, ft, 
                                                        strand, coord[0],
                                                        coord[1], score, ".", 
                                                        attributes)))
                                                elif n == 2:
                                                    (self.chrs[ch][parent]
                                                    .transcripts[pseudo_t]
                                                    .temp_UTRs.append(UTR
                                                            (ID, ch, source, ft, 
                                                            strand, coord[0],
                                                            coord[1], score, 
                                                            ".", attributes)))
                                                else:
                                                    (self.chrs[ch][parent]
                                                    .transcripts[pseudo_t].miRNAs
                                                    .append(Feature(ID, ch, 
                                                                    source, ft,
                                                                    strand,
                                                                    coord[0],
                                                                    coord[1],
                                                                    score,  ".",
                                                                    attributes)))                                                  
                                            else:
                                                print(f"{self.id} Error: {parent} "
                                                    "pseudogene already had a "
                                                    f"transcript which {ft} "
                                                    f"subfeature {ID} ignores")

                                        # subset of cases where subfeature directly
                                        # points to a gene that is not a pseudogene
                                        else:
                                            # gene without transcripts pointed to
                                            if (self.chrs[ch][parent]
                                                .transcripts) == {}:
                                                print(f"{self.id} Error: {ft} "
                                                    f"subfeature {ID} references"
                                                    f" {parent} gene which is "
                                                    "not a pseudogene and has "
                                                    "no transcripts")
                                                (self.errors["subfeature_to_gene"]
                                                .append(ID))
                                            # correctly linking the subfeature to 
                                            # the single transcript that exists
                                            elif (len(self.chrs[ch][parent]
                                                    .transcripts) == 1):
                                                temp_id = list(self.chrs[ch]
                                                            [parent].transcripts
                                                            .keys())[0]
                                                if n == 4:
                                                    (self.chrs[ch][parent]
                                                    .transcripts[temp_id]
                                                    .temp_CDSs.append(Feature
                                                            (ID, ch, source, ft, 
                                                            strand, coord[0],
                                                            coord[1], score,
                                                            phase, attributes)))
                                                elif n == 3:
                                                    (self.chrs[ch][parent]
                                                    .transcripts[temp_id].exons
                                                    .append
                                                    (Exon(ID, ch, source, ft, 
                                                        strand, coord[0],
                                                        coord[1], score, ".", 
                                                        attributes)))
                                                elif n == 2:
                                                    (self.chrs[ch][parent]
                                                    .transcripts[temp_id]
                                                    .temp_UTRs.append(UTR
                                                            (ID, ch, source, ft, 
                                                            strand, coord[0],
                                                            coord[1], score,
                                                            ".", attributes)))
                                                else:
                                                    (self.chrs[ch][parent]
                                                    .transcripts[temp_id]
                                                    .miRNAs.append(Feature
                                                            (ID, ch, source, ft, 
                                                            strand, coord[0],
                                                            coord[1], score,
                                                             ".", attributes)))
                                            else:
                                                print(f"{self.id} Error: {ft} "
                                                    f"subfeature {ID} references"
                                                    f"{parent} gene which is not"
                                                    "a pseudogene and has "
                                                    "multiple transcripts")
                                                (self.errors["subfeature_to_gene"]
                                                .append(ID))

                                    else:
                                        if self.liftoff:
                                            print(f"{self.id} Warning: {ft} "
                                                f"subfeature {ID} references "
                                                f"{parent} which is not found in"
                                                " the gff, possibly not "
                                                "transfered during liftoff")
                                            (self.warnings
                                            ["missing_subfeature_parent_liftoff"]
                                            .append(ID))
                                        else:
                                            print(f"{self.id} Warning: {ft} "
                                                f"subfeature {ID} references "
                                                f"{parent} which is not found in"
                                                " the gff")
                                            (self.warnings
                                            ["missing_subfeature_parent"]
                                            .append(ID))

                                elif created_genes:
                                    count = transcript_info[parent][3]
                                    if count > 1:
                                        for x in range(count):
                                            if f"{parent}_{x+1}_gene" in self.chrs[ch]:
                                                if f"{parent}_{x+1}" in self.chrs[ch][f"{parent}_{x+1}_gene"].transcripts:
                                                    start_parent = self.chrs[ch][f"{parent}_{x+1}_gene"].transcripts[f"{parent}_{x+1}"].start
                                                    end_parent = self.chrs[ch][f"{parent}_{x+1}_gene"].transcripts[f"{parent}_{x+1}"].end
                                                    if ((start_parent <= coord[0])
                                                        and (end_parent >= coord[1])):
                                                        if n == 4:
                                                            # CDS segments
                                                            (self.chrs[ch][f"{parent}_{x+1}_gene"].transcripts[f"{parent}_{x+1}"].temp_CDSs
                                                            .append(Feature(ID, ch, source,
                                                                            ft, strand, 
                                                                            coord[0], coord[1],
                                                                            score, phase, 
                                                                            f"ID={ID},Parent={ID}_{x+1}")))
                                                        
                                                        elif n == 3:
                                                            (self.chrs[ch][f"{parent}_{x+1}_gene"].transcripts[f"{parent}_{x+1}"].exons
                                                            .append(Exon(ID, ch, source, ft, 
                                                                        strand, coord[0], 
                                                                        coord[1], score, ".", 
                                                                        f"ID={ID},Parent={ID}_{x+1}")))
                                                        
                                                        elif n == 2:
                                                            (self.chrs[ch][f"{parent}_{x+1}_gene"].transcripts[f"{parent}_{x+1}"].temp_UTRs
                                                            .append(UTR(ID, ch, source, ft,
                                                                        strand, coord[0], 
                                                                        coord[1], score, ".", 
                                                                        f"ID={ID},Parent={ID}_{x+1}")))
                                                        else:
                                                            (self.chrs[ch][f"{parent}_{x+1}_gene"].transcripts[f"{parent}_{x+1}"].miRNAs
                                                            .append(Feature(ID, ch, source,
                                                                            ft, strand,
                                                                            coord[0], coord[1],
                                                                            score,  ".",
                                                                            f"ID={ID},Parent={ID}_{x+1}")))
                                    else:
                                        for g in self.chrs[ch].values():
                                            if parent in g.transcripts:
                                                if n == 4:
                                                    # CDS segments
                                                    (g.transcripts[parent].temp_CDSs
                                                    .append(Feature(ID, ch, source,
                                                                    ft, strand, 
                                                                    coord[0], coord[1],
                                                                    score, phase, 
                                                                    attributes)))
                                                
                                                elif n == 3:
                                                    (g.transcripts[parent].exons
                                                    .append(Exon(ID, ch, source, ft, 
                                                                strand, coord[0], 
                                                                coord[1], score, ".", 
                                                                attributes)))
                                                
                                                elif n == 2:
                                                    (g.transcripts[parent].temp_UTRs
                                                    .append(UTR(ID, ch, source, ft,
                                                                strand, coord[0], 
                                                                coord[1], score, ".", 
                                                                attributes)))
                                                else:
                                                    (g.transcripts[parent].miRNAs
                                                    .append(Feature(ID, ch, source,
                                                                    ft, strand,
                                                                    coord[0], coord[1],
                                                                    score,  ".",
                                                                    attributes)))                                      
                                else:
                                    found = 0
                                    # looking for parent amongst transcripts
                                    for g in self.chrs[ch].values():
                                        if parent in g.transcripts:
                                            found += 1
                                            if n == 4:
                                                # CDS segments
                                                (g.transcripts[parent].temp_CDSs
                                                .append(Feature(ID, ch, source,
                                                                ft, strand, 
                                                                coord[0], coord[1],
                                                                score, phase, 
                                                                attributes)))
                                            
                                            elif n == 3:
                                                (g.transcripts[parent].exons
                                                .append(Exon(ID, ch, source, ft, 
                                                            strand, coord[0], 
                                                            coord[1], score, ".", 
                                                            attributes)))
                                            
                                            elif n == 2:
                                                (g.transcripts[parent].temp_UTRs
                                                .append(UTR(ID, ch, source, ft,
                                                            strand, coord[0], 
                                                            coord[1], score, ".", 
                                                            attributes)))
                                            else:
                                                (g.transcripts[parent].miRNAs
                                                .append(Feature(ID, ch, source,
                                                                ft, strand,
                                                                coord[0], coord[1],
                                                                score,  ".",
                                                                attributes)))                                    
                                    if found > 1:
                                        print(f"{self.id} Error: same {parent} transcript from {ft} {ID} found for different genes.")
                                        self.errors["repeat_transcript_different_genes"].append(ID)          
                    else:
                        self.atypical_features.append((Feature(ID, ch, source, ft, strand, coord[0], coord[1], score, ".", attributes)))
                    progress_bar.update(1)
                progress_bar.close()

        now = time.time()
        lapse = now - start

        print(f"\nCreating {self.id} annotation object took {round(lapse/60, 1)} minutes\n")
        
        if self.features == ["nucleotide_to_protein_match"]:
            self.update(original_annotation=original_annotation, genome=genome, sort_processes=sort_processes, define_synteny=define_synteny)
        else:
            self.update(original_annotation=original_annotation, genome=genome, sort_processes=sort_processes, define_synteny=define_synteny)
            
    
        if "CDS" not in self.features and rework_CDSs:
            self.rework_CDSs(genome)
            self.update(original_annotation=original_annotation, genome=genome, sort_processes=sort_processes, define_synteny=define_synteny)

    def copy(self):
        return copy.deepcopy(self)
    
    def update(self, original_annotation:object=None, rename_ids_minimal=False, extra_attributes:bool=False, genome:object=None, define_synteny:bool=False, sort_processes:int=2):
        start = time.time()
        # Check if stdout or stderr are redirected to files
        stdout_redirected = not sys.stdout.isatty()
        stderr_redirected = not sys.stderr.isatty()

        # Disable tqdm if stdout or stderr are redirected
        if stdout_redirected or stderr_redirected:
            disable = True
        else:
            disable = False
        progress_bar = tqdm(total=len(self.all_gene_ids.keys()), disable=disable,
                                bar_format=(
                    f'\033[1;62mUpdating {self.id} genes:\033[0m '
                    '{percentage:3.0f}%|'
                    f'\033[1;62m{{bar}}\033[0m| '
                    '{n}/{total} [{elapsed}<{remaining}]'))
        for genes in self.chrs.values():
            for g in genes.values():
                progress_bar.update(1)
                for t in g.transcripts.values():
                    t.update()
                    if t.polycistronic == "maybe":
                        self.warnings["possible_policistronic_transcript"].append(t.id) 
                    elif t.polycistronic == "yes":
                        self.warnings["multiple_CDSs_per_transcript"].append(t.id)
                g.update()
        self.update_features()
        progress_bar.close()
        self.rename_ids(minimal=rename_ids_minimal, extra_attributes=extra_attributes)
        self.remove_missing_transcript_parent_references(extra_attributes=extra_attributes)
        self.correct_gene_transcript_and_subfeature_coordinates()
        if not self.sorted:
            self.sort_genes(processes=sort_processes)
        if define_synteny:
            self.define_synteny(original_annotation=original_annotation, sort_processes=sort_processes)
        if self.liftoff and original_annotation != None:
            for g_id in original_annotation.all_gene_ids:
                #improve at some point, slow with massive lists
                if g_id not in self.all_gene_ids:
                    self.unmapped.append(g_id)
        self.update_stats(genome=genome)

        self.update_suffixes()

        now = time.time()
        lapse = now - start
        print(f"\nWhole update process for {self.id} annotation object took {round(lapse/60, 1)} minutes\n")

    def update_suffixes(self):
        #print(f"\nUpdating suffixes for {self.id}")
        self.feature_suffix = ""
        self.suffix = ""

        if self.combined:
            self.feature_suffix += "_combined"
        if self.small_cds_removed:
            self.feature_suffix += "_minus_small_CDSs"
        if self.coding_removed:
            self.feature_suffix += "_minus_coding"
        if self.non_coding_removed:
            self.feature_suffix += "_minus_non_coding"
        if self.transposable_removed:
            self.feature_suffix += "_minus_TE"
        if self.non_transposable_removed:
            self.feature_suffix += "_minus_non_TE"

        self.suffix = self.feature_suffix

        if self.confrenamed:
            self.suffix += "_confrenamed"
        if self.aegis:
            self.suffix += "_aegis"
        if self.featurecounts:
            self.suffix += "_fcounts"
        if self.dapmod:
            self.suffix += "_dapmod"
        elif self.dapfit:
            self.suffix += "_dapfit"
        if self.symbols_added:
            self.suffix += "_plus_symbols"
        #print(f"\nUpdated suffixes for {self.id}")

    def update_features(self, standardise=True):
        #print(f"\nUpdating features for {self.id}")
        if standardise:
            for genes in self.chrs.values():
                for g in genes.values():
                    g.feature = Annotation.features_r[g.feature]
                    for t in g.transcripts.values():
                        if t.feature not in Annotation.noncoding_transcripts:
                            t.feature = "mRNA"
                        for e in t.exons:
                            e.feature = "exon"
        new_features = set()
        for genes in self.chrs.values():
            for g in genes.values():
                new_features.add(g.feature)
                for t in g.transcripts.values():
                    new_features.add(t.feature)
                    for c in t.CDSs.values():
                        new_features.add(c.feature)
                    for e in t.exons:
                        new_features.add(e.feature)
        self.features = list(new_features)
        self.features.sort()
        #print(f"\nUpdated features for {self.id}")

    def mark_transposable_element_genes(self, TE_genes_file):
        TE_genes = set()
        f_in = open(TE_genes_file)
        for line in f_in:
            line = line.strip()
            if line != "":
                TE_genes.add(line)
        f_in.close()
        for genes in self.chrs.values():
            for g in genes.values():
                if g.id in TE_genes:
                    g.transposable = True
        self.update()

    def mark_rRNA_transcripts(self, rRNA_transcripts_file, clean:bool=True):
        rRNA_transcripts = set()
        f_in = open(rRNA_transcripts_file)
        for line in f_in:
            line = line.strip()
            if line != "":
                rRNA_transcripts.add(line)
        f_in.close()
        for genes in self.chrs.values():
            for g in genes.values():
                for t in g.transcripts.values():
                    if t.id in rRNA_transcripts:
                        t.feature = "rRNA"
                        t.CDSs = {}
                        t.coding = False
                        t.frame = "."
                        t.blast_hits = []
                        t.main = False
                        t.temp_CDSs = []
                        t.temp_UTRs = []
        if clean:
            self.remove_other_mRNA_transcripts_from_rRNA_genes()
        self.update(rename_ids_minimal=True)
    
    def remove_other_mRNA_transcripts_from_rRNA_genes(self):
        for chrom, genes in self.chrs.items():
            for g in genes.values():
                rRNA_positive = False
                mRNA_positive = False
                for t in g.transcripts.values():
                    if t.feature == "rRNA":
                        rRNA_positive = True
                    elif t.feature == "mRNA":
                        mRNA_positive = True
                if rRNA_positive == True and mRNA_positive == True:
                    mRNA_transcripts_to_remove = []
                    for t in g.transcripts.values():
                        if t.feature == "mRNA":
                            mRNA_transcripts_to_remove.append(t.id)
                    for t_id in mRNA_transcripts_to_remove:
                        del self.chrs[chrom][g.id].transcripts[t_id]   

    def correct_gene_transcript_and_subfeature_coordinates(self):
        #print(f"Correcting feature coordinates for {self.id}")

        # fixing exon/CDS sizes
        for genes in self.chrs.values():
            for g in genes.values():
                for t in g.transcripts.values():
                    for c in t.CDSs.values():
                        if c.start < t.exons[0].start:
                            print(f"Warning: {c.id} start should not be earlier than for first {t.id} exon, proceeding to fix {self.id}")
                            t.exons[0].start = c.start
                            t.exons[0].update_size()
                            self.sorted = False
                        if c.end > t.exons[-1].end:
                            print(f"Warning: {c.id} end should not extend beyond the last {t.id} exon, proceeding to fix {self.id}")
                            t.exons[-1].end = c.end
                            t.exons[-1].update_size()
                            self.sorted = False

        # fixing transcript/exon sizes
        for genes in self.chrs.values():
            for g in genes.values():
                for t in g.transcripts.values():
                    if t.exons != []:
                        if t.exons[0].start < t.start:
                            print(f"First exon start should not be earlier than for {t.id}, proceeding to fix {self.id}")
                        elif t.exons[0].start > t.start:
                            print(f"First exon should not start later than {t.id}, proceeding to fix {self.id}")
                        if t.exons[-1].end < t.end:
                            print(f"Last exon should not finish earlier than {t.id}, proceeding to fix {self.id}")
                        elif t.exons[-1].end > t.end:
                            print(f"Last exon should not finish later than {t.id}, proceeding to fix {self.id}")
                        if t.start != t.exons[0].start or t.end != t.exons[-1].end:
                            t.start = t.exons[0].start
                            t.end = t.exons[-1].end
                            t.update_size()
                            self.sorted = False 

        # fixing gene/transcript sizes
        for genes in self.chrs.values():
            for g in genes.values():
                earliest_start = None
                latest_end = None
                for n, t in enumerate(g.transcripts.values()):
                    if t.start < g.start:
                        print(f"{t.id} start should not be earlier than for {g.id}, proceeding to fix {self.id}")
                        g.start = t.start
                        self.sorted = False
                    if t.end > g.end:
                        print(f"{t.id} end should not extend beyond {g.id}, proceeding to fix {self.id}")
                        g.end = t.end
                        self.sorted = False
                    if n == 0:
                        earliest_start = t.start
                        latest_end = t.end
                    else:
                        if t.start < earliest_start:
                            earliest_start = t.start
                        if t.end > latest_end:
                            latest_end = t.end
                if earliest_start != None:
                    if g.start != earliest_start or g.end != latest_end:
                        print(f"{g.id} was too long and had to be trimmed to longest transcript ({self.id})")
                        g.start = earliest_start
                        g.end = latest_end
                        g.update_size()
                        self.sorted = False
        #print(f"Corrected feature coordinates for {self.id}")

    def generate_sequences(self, genome:object, just_CDSs:bool=False):
        start = time.time()
        for o in self.atypical_features:
            o.generate_sequence(genome)
        if not just_CDSs:
            for genes in self.chrs.values():
                for g in genes.values():
                    g.generate_sequence(genome)
                    for t in g.transcripts.values():
                        t.generate_sequence(genome)
                        for c in t.CDSs.values():
                            c.generate_sequence(genome)
            self.generated_all_sequences = True
            self.contains_all_sequences = True 
        else:
            for genes in self.chrs.values():
                for g in genes.values():
                    for t in g.transcripts.values():
                        for c in t.CDSs.values():
                            c.generate_sequence(genome)
        self.generated_CDS_sequences = True
        self.contains_CDS_sequences = True
        self.generated_protein_sequences = True
        self.contains_protein_sequences = True

        now = time.time()
        lapse = now - start
        print(f"\nGenerating sequences for {self.id} annotation object "
              f"took {round(lapse/60, 1)} minutes\n")

    def clear_sequences(self, just_hard=False, keep_proteins:bool=False):
        start = time.time()
        for o in self.atypical_features:
            o.clear_sequence(just_hard=just_hard)
        for genes in self.chrs.values():
            for g in genes.values():
                g.clear_sequence(just_hard=just_hard)
                for t in g.transcripts.values():
                    t.clear_sequence(just_hard=just_hard)
                    for c in t.CDSs.values():
                        c.clear_sequence(just_hard=just_hard, keep_proteins=keep_proteins)
        if not keep_proteins:
            self.contains_protein_sequences = False
        self.contains_CDS_sequences = False
        self.contains_all_sequences = False

        now = time.time()
        lapse = now - start
        print(f"Clearing sequences of {self.id} annotation object took "
              f"{round(lapse, 1)} seconds\n")

    def generate_promoters(self, genome:object, promoter_size:int=2000, from_ATG:bool=False, count_prom_size_from_TSS:bool=True):
        """
        - from_ATG will only be paid attention if CDS exists.
        - count_prom_size_from_TSS can only be set to false if from_ATG is True
            - in this case the promoter will start counting from the ATG
                - ortherwise the upstream ATG region is added to std. promoter
        """
        self.promoter_types = "standard"
        self.promoter_size = promoter_size
        if from_ATG and count_prom_size_from_TSS:
            self.promoter_types = "plus_5_UTR"
        elif from_ATG and not count_prom_size_from_TSS:
            self.promoter_types = "just_upstream_ATG"

        for genes in self.chrs.values():
            for g in genes.values():
                for t in g.transcripts.values():
                    t.generate_promoter(promoter_size, from_ATG, count_prom_size_from_TSS)

        self.generate_sequences(genome)

    def find_motifs(self, query_genes:list, motif:str, motif_length, glistname, motifid, tfname, backlist:list=[], backlistname:str="", custom_path:str=""):
        # Check if stdout or stderr are redirected to files
        stdout_redirected = not sys.stdout.isatty()
        stderr_redirected = not sys.stderr.isatty()

        # Disable tqdm if stdout or stderr are redirected
        if stdout_redirected or stderr_redirected:
            disable = True
        else:
            disable = False

        bin_division = 30
        bins_genome_division = 30
        random_ids = self.return_random_gene_ids(len(query_genes), to_avoid=query_genes)
        if backlist == []:
            total = (len(query_genes) * 2) + len(self.all_gene_ids.keys())
        else:
            total = (len(query_genes) * 2) + len(backlist)
        progress_bar = tqdm(total=total, disable=disable,
                                bar_format=(
                    f'\033[1;94;1mScanning {glistname} genes for {tfname} {motif} ({motifid}):\033[0m '
                    '{percentage:3.0f}%|'
                    f'\033[1;94;1m{{bar}}\033[0m| '
                    '{n}/{total} [{elapsed}<{remaining}]'))
        if custom_path == "":
            output_path = self.path + "tf_motifs/"
            output_file = output_path 
        else:
            output_file = custom_path
        if output_file[-1] != "/":
            output_file += "/"
        common_multifasta = self.path + "tf_motifs/"
        system(f"mkdir -p {common_multifasta}")
        output_file += f"{tfname}/"
        system(f"mkdir -p {output_file}")

        output_file += f"{tfname}_{motifid}"

        motif = motif.upper()
        promoter_length = 0

        for genes in self.chrs.values():
            for g in genes.values():
                for t in g.transcripts.values():
                    if t.main:
                        if hasattr(t, "promoter"):
                            if t.promoter.seq != "":
                                promoter_length = len(t.promoter.seq)
                                break
                break
            break
                    
        multifasta_out = ""
        # critical thing to understand here is that we are looking at how a motif
        # is oriented with regards to the TSS
        towards_occurrences = []
        against_occurrences = []
        interest_proportion = 0
        avg_motifs_interest = []
        for id in query_genes:
            progress_bar.update(1)
            ch = self.all_gene_ids[id]
            g = self.chrs[ch][id]
            for t in g.transcripts.values():
                if t.main:
                    p = t.promoter
                    multifasta_out += f">{g.id}_promoter\n{p.seq}\n"
                    occurrences_t = find_all_occurrences(motif, p.seq)
                    occurrences_a = find_all_occurrences(motif, reverse_complement(p.seq))
                    occurrence_count_total = len(occurrences_t) + len(occurrences_a)
                    if occurrence_count_total != 0:
                        avg_motifs_interest.append(occurrence_count_total)
                    if occurrences_t != [] or occurrences_a != []:
                        interest_proportion += 1
                    for m in occurrences_t:
                        towards_occurrences.append((m[0], m[1], p.start + m[0], p.start + m[1], m[2], "same", g.id, g.name, g.strand, g.ch))
                    for m in occurrences_a:
                        against_occurrences.append((m[0], m[1], p.end - m[0], p.end - m[1], m[2], "different", g.id, g.name, g.strand, g.ch))


        midpoints = []
        for o in towards_occurrences:
            midpoint = (o[0] + o[1]) // 2
            midpoint -= promoter_length
            midpoints.append(midpoint)

        for o in against_occurrences:
            midpoint = (o[0] + o[1]) // 2
            midpoint += 1
            midpoint = (midpoint * -1)
            midpoints.append(midpoint)
        
        if backlist == []:
            f_out = open(f"{common_multifasta}{glistname}_{promoter_length/1000}kb.fasta", "w", encoding="utf-8")
            f_out.write(multifasta_out)
            f_out.close()

        all_occurrences = towards_occurrences + against_occurrences
        df = pd.DataFrame(all_occurrences, columns=["start", "end", "genomic_start", "genomic_end", "sequence", "orientation with respect to gene", "gene_id", "gene_name", "gene_strand", "chromosome"])
        if backlist == []:
            df.to_csv(f"{output_file}_{glistname}.csv", sep="\t")
        
        interest_count = len(midpoints)
        plt.hist(midpoints, bins=(promoter_length//bin_division), color='skyblue', edgecolor='skyblue')
        plt.grid(b=False)
        plt.title(f"{self.id} {tfname} {motifid} {glistname} histogram\npromoters with motif: ({interest_proportion}/{len(query_genes)})")
        plt.xlabel("promoter position")
        plt.ylabel(f"motif occurrence count (total: {interest_count})")
        plt.grid(True)
        if backlist == []:
            plt.savefig(f"{output_file}_{glistname}.pdf")
        plt.close()

        multifasta_out = ""
        # critical thing to understand here is that we are looking at how a motif
        # is oriented with regards to the TSS
        towards_occurrences = []
        against_occurrences = []
        avg_motifs_random = []
        random_proportion = 0
        for id in random_ids:
            progress_bar.update(1)
            ch = self.all_gene_ids[id]
            g = self.chrs[ch][id]
            for t in g.transcripts.values():
                if t.main:
                    p = t.promoter
                    multifasta_out += f">{g.id}_promoter\n{p.seq}\n"
                    occurrences_t = find_all_occurrences(motif, p.seq)
                    occurrences_a = find_all_occurrences(motif, reverse_complement(p.seq))
                    if occurrence_count_total != 0:
                        avg_motifs_random.append(occurrence_count_total)
                    if occurrences_t != [] or occurrences_a != []:
                        random_proportion += 1
                    for m in occurrences_t:
                        towards_occurrences.append((m[0], m[1], p.start + m[0], p.start + m[1], m[2], "same", g.id, g.name, g.strand, g.ch))
                    for m in occurrences_a:
                        against_occurrences.append((m[0], m[1], p.end - m[0], p.end - m[1], m[2], "different", g.id, g.name, g.strand, g.ch))

        midpoints = []
        for o in towards_occurrences:
            midpoint = (o[0] + o[1]) // 2
            midpoint -= promoter_length
            midpoints.append(midpoint)

        for o in against_occurrences:
            midpoint = (o[0] + o[1]) // 2
            midpoint += 1
            midpoint = (midpoint * -1)
            midpoints.append(midpoint)
        
        if backlist == []:
            f_out = open(f"{output_file}_{glistname}_random.fasta", "w", encoding="utf-8")
            f_out.write(multifasta_out)
            f_out.close()

        all_occurrences = towards_occurrences + against_occurrences
        df = pd.DataFrame(all_occurrences, columns=["start", "end", "genomic_start", "genomic_end", "sequence", "orientation with respect to gene", "gene_id", "gene_name", "gene_strand", "chromosome"])
        if backlist == []:
            df.to_csv(f"{output_file}_{glistname}_random.csv", sep="\t")

        random_count = len(midpoints)
        plt.hist(midpoints, bins=(promoter_length//bin_division), color='grey', edgecolor='grey')
        plt.grid(b=False)
        plt.title(f"{self.id} {tfname} {motifid} {glistname} random histogram\npromoters with motif: ({random_proportion}/{len(query_genes)})")
        plt.xlabel("promoter position")
        plt.ylabel(f"motif occurrence count (total: {random_count})")
        plt.grid(True)
        if backlist == []:
            plt.savefig(f"{output_file}_{glistname}_random.pdf")
        plt.close()

        # critical thing to understand here is that we are looking at how a motif
        # is oriented with regards to the TSS
        towards_occurrences = []
        against_occurrences = []
        avg_motifs_genomic = []
        genomic_proportion = 0
        if backlist == []:
            for id in self.all_gene_ids:
                progress_bar.update(1)
                ch = self.all_gene_ids[id]
                g = self.chrs[ch][id]
                for t in g.transcripts.values():
                    if t.main:
                        p = t.promoter
                        occurrences_t = find_all_occurrences(motif, p.seq)
                        occurrences_a = find_all_occurrences(motif, reverse_complement(p.seq))
                        occurrence_count_total = len(occurrences_t) + len(occurrences_a)
                        if occurrence_count_total != 0:
                            avg_motifs_genomic.append(occurrence_count_total)
                        if occurrences_t != [] or occurrences_a != []:
                            genomic_proportion += 1
                        for m in occurrences_t:
                            towards_occurrences.append((m[0], m[1], p.start + m[0], p.start + m[1], m[2], "same", g.id, g.name, g.strand, g.ch))
                        for m in occurrences_a:
                            against_occurrences.append((m[0], m[1], p.end - m[0], p.end - m[1], m[2], "different", g.id, g.name, g.strand, g.ch))

            midpoints = []
            for o in towards_occurrences:
                midpoint = (o[0] + o[1]) // 2
                midpoint -= promoter_length
                midpoints.append(midpoint)

            for o in against_occurrences:
                midpoint = (o[0] + o[1]) // 2
                midpoint += 1
                midpoint = (midpoint * -1)
                midpoints.append(midpoint)

            genomic_count = len(midpoints)
            plt.hist(midpoints, bins=(promoter_length//bins_genome_division), color='grey', edgecolor='grey')
            plt.grid(b=False)
            plt.title(f"{self.id} {tfname} {motifid} full genome histogram\npromoters with motif: ({genomic_proportion}/{len(self.all_gene_ids.keys())})")
            plt.xlabel("promoter position")
            plt.ylabel(f"motif occurrence count (total: {genomic_count})")
            plt.grid(True)
            plt.savefig(f"{output_file}_whole_genome.pdf")
            plt.close()
        
        else:
            for id in backlist:
                progress_bar.update(1)
                ch = self.all_gene_ids[id]
                g = self.chrs[ch][id]
                for t in g.transcripts.values():
                    if t.main:
                        p = t.promoter
                        occurrences_t = find_all_occurrences(motif, p.seq)
                        occurrences_a = find_all_occurrences(motif, reverse_complement(p.seq))
                        occurrence_count_total = len(occurrences_t) + len(occurrences_a)
                        if occurrence_count_total != 0:
                            avg_motifs_genomic.append(occurrence_count_total)
                        if occurrences_t != [] or occurrences_a != []:
                            genomic_proportion += 1
                        for m in occurrences_t:
                            towards_occurrences.append((m[0], m[1], p.start + m[0], p.start + m[1], m[2], "same", g.id, g.name, g.strand, g.ch))
                        for m in occurrences_a:
                            against_occurrences.append((m[0], m[1], p.end - m[0], p.end - m[1], m[2], "different", g.id, g.name, g.strand, g.ch))

            midpoints = []
            for o in towards_occurrences:
                midpoint = (o[0] + o[1]) // 2
                midpoint -= promoter_length
                midpoints.append(midpoint)

            for o in against_occurrences:
                midpoint = (o[0] + o[1]) // 2
                midpoint += 1
                midpoint = (midpoint * -1)
                midpoints.append(midpoint)

            genomic_count = len(midpoints)
            plt.hist(midpoints, bins=(promoter_length//bins_genome_division), color='grey', edgecolor='grey')
            plt.grid(b=False)
            plt.title(f"{self.id} {tfname} {motifid} {backlistname} as background histogram\npromoters with motif: ({genomic_proportion}/{len(backlist)})")
            plt.xlabel("promoter position")
            plt.ylabel(f"motif occurrence count (total: {genomic_count})")
            plt.grid(True)
            plt.savefig(f"{output_file}_{backlistname}_as_background.pdf")
            plt.close()

        progress_bar.close()

        # Counts of non-motif occurrences
        interest_non_count = (len(query_genes) * (int(promoter_length/motif_length)) * 2) - interest_count
        if backlist == []:
            genomic_non_count = (len(self.all_gene_ids.keys()) * (int(promoter_length/motif_length)) * 2) - genomic_count
        else:
            genomic_non_count = (len(backlist) * (int(promoter_length/motif_length)) * 2) - genomic_count
        interest_non_proportion = len(query_genes) - interest_proportion
        if backlist == []:
            genomic_non_proportion = len(self.all_gene_ids.keys()) - genomic_proportion
        else:
            genomic_non_proportion = len(backlist) - genomic_proportion

        odds_ratio_occurrences, p_value_occurrences = fisher_exact([[interest_count, genomic_count], 
                                                                    [interest_non_count, genomic_non_count]])

        odds_ratio_proportion, p_value_proportion = fisher_exact([[interest_proportion, genomic_proportion], 
                                                                  [interest_non_proportion, genomic_non_proportion]])

        promoter_percentage_interest = (interest_proportion / len(query_genes)) * 100
        if backlist == []:
            promoter_percentage_genome = (genomic_proportion / len(self.all_gene_ids.keys())) * 100
        else:
            promoter_percentage_genome = (genomic_proportion / len(backlist)) * 100

        return interest_count, genomic_count, p_value_occurrences, odds_ratio_occurrences, promoter_percentage_interest, promoter_percentage_genome, p_value_proportion, odds_ratio_proportion, avg_motifs_interest, avg_motifs_genomic, output_file

    def return_random_gene_ids(self, number:int=1, to_avoid:list=[], coding:bool=True):
        random_ids = []
        while len(random_ids) < number:
            r = random.choice(list(self.all_gene_ids.keys()))
            ch = self.all_gene_ids[r]
            if coding:
                if self.chrs[ch][r].coding and r not in random_ids and r not in to_avoid:
                    random_ids.append(r)
            elif r not in random_ids and r not in to_avoid:
                random_ids.append(r)

        return random_ids

    def export_all_features(self, feature_output:str="main", promoters:bool=True, verbose:bool=True, path:str="", most_specific_id_level="promoter"):
        """
        The "output" parameter can be both, main or all. This parameter only 
        affects promoter, transcript, CDS and protein sequences. If both is selected
        a "main features" file and and "all features" file will be produced.
        """

        feature_output_choices = ["main", "all", "both"]
        if feature_output not in feature_output_choices:
            raise ValueError(f"feature_output={feature_output} is not amongst the feature_output_choices={feature_output_choices} to export all features.")

        most_specific_id_level_choices = ["gene", "transcript", "CDS", "protein", "promoter"]

        if most_specific_id_level not in most_specific_id_level_choices:
            raise ValueError(f"most_specific_id_level={most_specific_id_level} is not amongst the most_specific_id_level_choices={most_specific_id_level_choices} to export proteins.")

        start = time.time()

        self.export_genes(verbose, path)

        if feature_output == "both":
            modes = [True, False]
        elif feature_output == "main":
            modes = [True]
        elif feature_output == "all":
            modes = [False]

        for b in modes:

            if most_specific_id_level == "promoter" or most_specific_id_level == "protein":
                self.export_proteins(b, verbose, path)
                self.export_CDSs(b, verbose, path)
                self.export_transcripts(b, verbose, path)
                if promoters:
                    self.export_promoters(b, verbose, path)

            elif most_specific_id_level == "CDS":
                self.export_proteins(b, verbose, path, used_id="CDS")
                self.export_CDSs(b, verbose, path)
                self.export_transcripts(b, verbose, path)
                if promoters:
                    self.export_promoters(b, verbose, path)

            elif most_specific_id_level == "transcript":
                self.export_proteins(b, verbose, path, used_id="transcript")
                self.export_CDSs(b, verbose, path, used_id="transcript")
                self.export_transcripts(b, verbose, path)
                if promoters:
                    self.export_promoters(b, verbose, path, used_id="transcript")

            elif most_specific_id_level == "gene":
                self.export_proteins(b, verbose, path, used_id="gene")
                self.export_CDSs(b, verbose, path, used_id="gene")
                self.export_transcripts(b, verbose, path, used_id="gene")
                if promoters:
                    self.export_promoters(b, verbose, path, used_id="gene")

        now = time.time()
        lapse = now - start
        print(f"Extracting {self.id} annotation features took "
              f"{round(lapse, 1)} seconds\n")
            
    def export_proteins(self, only_main:bool=True, verbose:bool=True, custom_path:str="", used_id:str="protein", unique_proteins_per_gene:bool=False):
        """
        Main proteins means only proteins obtained from the main CDSs of the
        main transcripts. This equates to one protein per gene.

        Verbose will include protein tag details and readthrough options, strand, chromosome, and coordinates.

        Proteins will be exported into fastas with their protein ids unless 
        CDS, transcript or gene has been selected as used_id. To avoid duplicate
        equal entries, when choosing gene or transcript only main proteins will
        be able to be output.

        valid_id_choices = ["gene", "transcript", "CDS", "protein"]
        """
        if custom_path == "":
            output_path = self.path + "features/"
            output_file = output_path 
        else:
            output_file = custom_path
        if output_file[-1] != "/":
            output_file += "/"
        system(f"mkdir -p {output_file}")
        output_file += self.id
        output_file += self.feature_suffix
        output_file += "_proteins"

        out = ""

        only_cds_main = only_main

        if used_id == "gene":
            only_main = True
            only_cds_main = True
            output_file += "_g_id_main"
        else:
            if used_id == "transcript":
                only_cds_main = True
                output_file += "_t_id"
            elif used_id == "CDS":
                output_file += "_c_id"
            elif used_id == "protein":
                output_file += "_p_id"

            if only_main:
                output_file += "_main"
            else:
                output_file += "_all"
                if unique_proteins_per_gene:
                    output_file += "_unique_per_gene"

        if verbose:
            output_file += "_coordinates"
            
        output_file += ".fasta"
        
        valid_id_choices = ["gene", "transcript", "CDS", "protein"]
        if used_id not in valid_id_choices:
            raise ValueError(f"used_id={used_id} is not amongst the valid_id_choices={valid_id_choices} to export proteins.")
        
        for genes in self.chrs.values():
            for g in genes.values():
                temp_cs = []
                for t in g.transcripts.values():
                    if only_main:
                        if t.main:
                            for c in t.CDSs.values():
                                if c.seq != "":
                                    if only_cds_main:
                                        if c.main:
                                            temp_cs.append(c)
                                    else:
                                        temp_cs.append(c)
                    else:
                        for c in t.CDSs.values():
                            if c.seq != "":
                                if only_cds_main:
                                    if c.main:
                                        temp_cs.append(c)
                                else:
                                    temp_cs.append(c)

                if not unique_proteins_per_gene:
                    final_cs = temp_cs.copy()

                else:
                    final_cs = []
                    if len(temp_cs) > 0:
                        final_cs.append(temp_cs[0])

                    for i, c1 in enumerate(temp_cs):
                        if i > 0:
                            add = True
                            for c2 in final_cs:
                                if c1.equal_segments(c2):
                                    add = False
                            if add:
                                final_cs.append(c1)

                for c in final_cs:

                    if used_id == "protein":
                        out += f">{c.protein.id}"
                    elif used_id == "CDS":
                        out += f">{c.id}"
                    elif used_id == "transcript":
                        out += f">{t.id}"
                    elif used_id == "gene":
                        out += f">{g.id}"

                    if c.protein.summary_tag and verbose:
                        out += f"|{c.protein.summary_tag}"
                    elif verbose:
                        out += f"|readthrough:{c.protein.readthrough}|{c.strand}|{c.ch}|{c.start}:{c.end}"

                    out += f"\n{c.protein.seq}\n"        
                     
        if out != "":
            f_out = open(output_file, "w", encoding="utf-8")
            f_out.write(out)
            f_out.close()
        else:
            print(f"Warning: Run self.generate_sequences(genome) on {self.id}")

    def export_unique_proteins(self, genome:object=None, custom_path:str=""):
        start = time.time()
        # Check if stdout or stderr are redirected to files
        stdout_redirected = not sys.stdout.isatty()
        stderr_redirected = not sys.stderr.isatty()

        # Disable tqdm if stdout or stderr are redirected
        if stdout_redirected or stderr_redirected:
            disable = True
        else:
            disable = False
        
        if custom_path == "":
            output_path = self.path + "features/"
            output_file = output_path 
        else:
            output_file = custom_path
        if output_file[-1] != "/":
            output_file += "/"

        if not self.contains_protein_sequences:
            if genome == None:
                print(f"You forgot to provide the genome for {self.id} and no sequences exist for CDSs")
            else:
                self.generate_sequences(genome)
                self.export_unique_proteins(custom_path=custom_path)
        else:
            system(f"mkdir -p {output_file}")
            output_file += f"{self.id}{self.feature_suffix}_unique_proteins.fasta"
            out = ""
            all_protein_seqs = {}
            self.all_protein_ids = {}
            for chrom, genes in self.chrs.items():
                for g in genes.values():
                    if g.coding:
                        for t in g.transcripts.values():
                            for c in t.CDSs.values():
                                if c.protein.seq != "":
                                    all_protein_seqs[c.protein.id] = c.protein.seq
                                    self.all_protein_ids[c.protein.id] = (chrom, g.id, t.id, c.id)

            progress_bar = tqdm(total=len(all_protein_seqs.keys()), disable=disable,
                            bar_format=(
                f'\033[1;91mDetermining and exporting unique {self.id} proteins:\033[0m '
                '{percentage:3.0f}%|'
                f'\033[1;91m{{bar}}\033[0m| '
                '{n}/{total} [{elapsed}<{remaining}]'))
            
            
            unique_sequences = {}
            self.protein_equivalences = {}

            for protein_id, sequence in all_protein_seqs.items():
                progress_bar.update(1)
                if sequence not in unique_sequences:
                    unique_sequences[sequence] = protein_id
                    self.protein_equivalences[protein_id] = []
                else:
                    first_protein_id = unique_sequences[sequence]
                    self.protein_equivalences[first_protein_id].append(protein_id)

            for sequence, protein_id in unique_sequences.items():
                out += f">{protein_id}\n{sequence}\n"

            f_out = open(output_file, "w", encoding="utf-8")
            f_out.write(out)
            f_out.close()

            progress_bar.close()

            now = time.time()
            lapse = now - start
            print(f"\nExporting unique {self.id} proteins took {round(lapse/60, 1)} minutes")

    def export_CDSs(self, only_main:bool=True, verbose:bool=True, custom_path:str="", used_id:str="CDS", unique_CDSs_per_gene:bool=False):
        """
        Main CDSs means only CDS sequence obtained from the main CDS of the
        main transcripts.

        Verbose will include strand, chromosome, and coordinates.

        CDSs will be exported into fastas with their CDS ids unless 
        transcript or gene has been selected as used_id. To avoid duplicate
        equal entries, when choosing gene or transcript only main CDSs will
        be able to be output.

        valid_id_choices = ["gene", "transcript", "CDS"]
        """
        if custom_path == "":
            output_path = self.path + "features/"
            output_file = output_path 
        else:
            output_file = custom_path
        if output_file[-1] != "/":
            output_file += "/"
        system(f"mkdir -p {output_file}")
        output_file += self.id
        output_file += self.feature_suffix
        output_file += "_CDSs"

        out = ""

        only_cds_main = only_main

        if used_id == "gene":
            only_main = True
            only_cds_main = True
            output_file += "_g_id_main"
        else:
            if used_id == "transcript":
                only_cds_main = True
                output_file += "_t_id"
            elif used_id == "CDS":
                output_file += "_c_id"

            if only_main:
                output_file += "_main"
            else:
                output_file += "_all"
                if unique_CDSs_per_gene:
                    output_file += "_unique_per_gene"

        if verbose:
            output_file += "_coordinates"
            
        output_file += ".fasta"
        
        valid_id_choices = ["gene", "transcript", "CDS"]
        if used_id not in valid_id_choices:
            raise ValueError(f"used_id={used_id} is not amongst the valid_id_choices={valid_id_choices} to export CDSs.")

        for genes in self.chrs.values():
            for g in genes.values():
                temp_cs = []
                for t in g.transcripts.values():
                    if only_main:
                        if t.main:
                            for c in t.CDSs.values():
                                if c.seq != "":
                                    if only_cds_main:
                                        if c.main:
                                            temp_cs.append(c)
                                    else:
                                        temp_cs.append(c)
                    else:
                        for c in t.CDSs.values():
                            if c.seq != "":
                                if only_main:
                                    if c.main:
                                        temp_cs.append(c)
                                else:
                                    temp_cs.append(c)

                if not unique_CDSs_per_gene:
                    final_cs = temp_cs.copy()
                else:
                    final_cs = []

                    if len(temp_cs) > 0:
                        final_cs.append(temp_cs[0])

                    for i, c1 in enumerate(temp_cs):
                        if i > 0:
                            add = True
                            for c2 in final_cs:
                                if c1.equal_segments(c2):
                                    add = False
                            if add:
                                final_cs.append(c1)

                for c in final_cs:

                    if used_id == "CDS":
                        out += f">{c.id}"
                    elif used_id == "transcript":
                        out += f">{t.id}"
                    elif used_id == "gene":
                        out += f">{g.id}"

                    if verbose:
                        out += f"|{c.strand}|{c.ch}|{c.start}:{c.end}"

                    out += f"\n{c.seq}\n"
                     
        if out != "":
            f_out = open(output_file, "w", encoding="utf-8")
            f_out.write(out)
            f_out.close()
        else:
            print(f"Warning: Run self.generate_sequences(genome) on {self.id}")

    def export_transcripts(self, only_main:bool=True, verbose:bool=True, custom_path:str="", used_id:str="transcript"):
        """
        Main means only main transcript sequences are exported.

        Verbose will include strand, chromosome, and coordinates.

        Transcripts will be exported into fastas with their transcript ids unless 
        gene has been selected as used_id. To avoid duplicate
        equal entries, when choosing gene only main transcripts will
        be able to be output.

        valid_id_choices = ["gene", "transcript"]
        """
        if custom_path == "":
            output_path = self.path + "features/"
            output_file = output_path 
        else:
            output_file = custom_path
        if output_file[-1] != "/":
            output_file += "/"
        system(f"mkdir -p {output_file}")
        output_file += self.id
        output_file += self.feature_suffix
        output_file += "_transcripts"

        out = ""

        if used_id == "gene":
            only_main = True
            output_file += "_g_id_main"
        elif used_id == "transcript":
            output_file += "_t_id"
            if only_main:
                output_file += "_main"
            else:
                output_file += "_all"

        if verbose:
            output_file += "_coordinates"
            
        output_file += ".fasta"

        valid_id_choices = ["gene", "transcript"]
        if used_id not in valid_id_choices:
            raise ValueError(f"used_id={used_id} is not amongst the valid_id_choices={valid_id_choices} to export transcripts.")
        
        for genes in self.chrs.values():
            for g in genes.values():
                for t in g.transcripts.values():
                    if t.seq != "":
                        if only_main and not t.main:
                            continue

                        if used_id == "transcript":    
                            out += f">{t.id}"
                        elif used_id == "gene":
                            out += f">{g.id}"

                        if verbose:
                            out += f"|{t.strand}|{t.ch}|{t.start}:{t.end}"
                        out += f"\n{t.seq}\n"

        if out != "":
            f_out = open(output_file, "w", encoding="utf-8")
            f_out.write(out)
            f_out.close()
        else:
            print(f"Warning: Run self.generate_sequences(genome) on {self.id}")

    def export_genes(self, verbose:bool=True, custom_path:str=""):
        if custom_path == "":
            output_path = self.path + "features/"
            output_file = output_path 
        else:
            output_file = custom_path
        if output_file[-1] != "/":
            output_file += "/"
        system(f"mkdir -p {output_file}")
        output_file += self.id
        output_file += self.feature_suffix
        output_file += "_genes"
        if verbose:
            output_file += "_coordinates"
        output_file += ".fasta"
        out = ""

        for genes in self.chrs.values():
            for g in genes.values():
                if g.seq != "":
                    out += f">{g.id}"
                    if verbose:
                        out += f"|{g.strand}|{g.ch}|{g.start}:{g.end}"
                    out += f"\n{g.seq}\n"
        if out != "":
            f_out = open(output_file, "w", encoding="utf-8")
            f_out.write(out)
            f_out.close()
        else:
            print(f"Warning: Run self.generate_sequences(genome) on {self.id}")

    def export_promoters(self, only_main:bool=True, verbose:bool=True, custom_path:str="", used_id:str="promoter"):
        """

        Verbose will include promoter type, strand, chromosome, and coordinates.

        """
        if custom_path == "":
            output_path = self.path + "features/"
            output_file = output_path 
        else:
            output_file = custom_path
        system(f"mkdir -p {output_file}")
        if output_file[-1] != "/":
            output_file += "/"

        output_file += self.id
        output_file += self.feature_suffix
        output_file += "_promoters"

        out = ""

        if used_id == "gene":
            only_main = True
            output_file += "_g_id_main"
        else:
            if used_id == "transcript":
                output_file += "_t_id"
            if used_id == "promoter":
                output_file += "_prom_id"
            if only_main:
                output_file += "_main"
            else:
                output_file += "_all"

        if verbose:
            output_file += "_coordinates"
            
        output_file += f"_{self.promoter_size}_{self.promoter_types}.fasta"

        valid_id_choices = ["gene", "transcript", "promoter"]

        if used_id not in valid_id_choices:
            raise ValueError(f"used_id={used_id} is not amongst the valid_id_choices={valid_id_choices} to export promoters.")

        for genes in self.chrs.values():
            for g in genes.values():
                for t in g.transcripts.values():
                    if only_main and not t.main:
                        continue

                    if hasattr(t, "promoter"):
                        if t.promoter.seq != "":

                            if used_id == "promoter":
                                out += f">{t.promoter.id}"
                            elif used_id == "transcript":
                                out += f">{t.id}"
                            elif used_id == "gene":
                                out += f">{g.id}"

                            if verbose:
                                out += f"|{t.promoter.type}|{t.strand}|{t.ch}|{t.promoter.start}:{t.promoter.end}"

                            out += f"\n{t.promoter.seq}\n"
                        else:
                            print(f"Warning: transcript {t.id} from annotation {self.id} has no promoter sequence.\n")
                    else:
                        print(f"Warning: transcript {t.id} from annotation {self.id} has no promoter.\n")

        if out != "":
            f_out = open(output_file, "w", encoding="utf-8")
            f_out.write(out)
            f_out.close()
        else:
            print(f"Warning: Run self.generate_promoters(genome) and self.generate_sequences(genome) on {self.id}")


    def combine_transcripts(self, genome:object, respect_non_coding:bool=False):
        for genes in self.chrs.values():
            for g in genes.values():
                g.combine_transcripts(genome, respect_non_coding=respect_non_coding)
        self.sorted = False
        self.update()
        self.combined = True
        self.update_suffixes()

    def sort_genes(self, processes:int=2):
        print(f"\nSorting genes for {self.id}")
        # Check if stdout or stderr are redirected to files
        stdout_redirected = not sys.stdout.isatty()
        stderr_redirected = not sys.stderr.isatty()

        # Disable tqdm if stdout or stderr are redirected
        if stdout_redirected or stderr_redirected:
            disable = True
        else:
            disable = False
        progress_bar = tqdm(total=len(self.all_gene_ids.keys()), disable=disable,
                        bar_format=(
            f'\033[38;2;210;180;140m\033[1mSorting {self.id} genes:\033[0m '
            '{percentage:3.0f}%|'
            f'\033[38;2;210;180;140m{{bar}}\033[0m| '
            '{n}/{total} [{elapsed}<{remaining}]'))

        if processes > 1:
            print(f"Parallel mode chosen for {self.id}")
            with Pool(processes=processes) as pool:
                print(f"Parallel sorting started pool for {self.id}")
                for chrom, sorted_genes in pool.starmap(sort_and_update_genes, [(chrom, genes) for chrom, genes in self.chrs.items()]):
                    #print(f"Parallel sorting {chrom} genes for {self.id}")
                    self.chrs[chrom] = sorted_genes
                    progress_bar.update(len(sorted_genes))
        else:
            print(f"Sequential mode chosen for {self.id}")
            for chrom, genes in self.chrs.items():
                sorted_genes = sorted(genes.values())
                #print(f"Sequential sorting {chrom} genes for {self.id}")
                self.chrs[chrom] = {g.id: g.copy() for g in sorted_genes}
                progress_bar.update(len(sorted_genes))

        progress_bar.close()
        self.sorted = True
        print(f"Sorted genes for {self.id}")

    def define_synteny(self, original_annotation:object, sort_processes:int=2):
        print(f"\nDefining synteny for {self.id} annotation genes")
        start = time.time()
        if not self.sorted:
            self.sort_genes(processes=sort_processes)
        # defining synteny
        for genes in self.chrs.values():
            gene_list = list(genes.values())
            for n, g in enumerate(gene_list):
                g.synteny_order = n
                # This works even if only a single gene has been annotated
                # in a chromosome or scaffold
                if len(genes) == 1:
                    g.previous_gene = False
                    g.next_gene = False
                elif n == 0:
                    g.previous_gene = False
                    g.next_gene = gene_list[n+1].id
                elif n != (len(genes) - 1):
                    g.previous_gene = gene_list[n-1].id
                    g.next_gene = gene_list[n+1].id
                else:
                    g.previous_gene = gene_list[n-1].id
                    g.next_gene = False

        if self.liftoff:
            for genes in self.chrs.values():
                for g_id, g in genes.items():
                    # this extra bit is for extra liftoff copies
                    if g_id not in original_annotation.all_gene_ids:
                        g_id = "_".join(g_id.split("_")[:-1])
                        if g_id not in original_annotation.all_gene_ids:
                            continue

                    g.old_previous_gene = original_annotation.chrs[original_annotation.all_gene_ids[g_id]][g_id].previous_gene
                    g.old_next_gene = original_annotation.chrs[original_annotation.all_gene_ids[g_id]][g_id].next_gene                
                    g.old_synteny_order = original_annotation.chrs[original_annotation.all_gene_ids[g_id]][g_id].synteny_order

                    if g.old_previous_gene == g.previous_gene and g.old_next_gene == g.next_gene:
                        g.conserved_synteny = True
                    else: 
                        g.conserved_synteny = False
        now = time.time()
        lapse = now - start
        print(f"\nDefining synteny for {self.id} annotation genes took {round(lapse, 1)} seconds\n")

    def calculate_transcript_masking(self, hard_masked_genome:object):
        for genes in self.chrs.values():
            for g in genes.values():
                g.generate_hard_sequence(hard_masked_genome)
                for t in g.transcripts.values():
                    t.generate_hard_sequence(hard_masked_genome)
                    t.calculate_masking()
                    t.clear_sequence(just_hard=True)
                    for c in t.CDSs.values():
                        c.generate_hard_sequence(hard_masked_genome)
                        c.calculate_masking()
                        c.clear_sequence(just_hard=True)
        self.update()
    
    def calculate_gc_content(self):
        for genes in self.chrs.values():
            for g in genes.values():
                g.calculate_gc_content()
                for t in g.transcripts.values():
                    t.calculate_gc_content()
                    for c in t.CDSs.values():
                        c.calculate_gc_content()

    def update_stats(self, custom_path:str="", export:bool=False, genome:object=None, max_x:int=None):
        #print(f"\nUpdating stats for {self.id}")
        if not self.generated_all_sequences or not self.contains_protein_sequences:
            if genome != None:
                self.generate_sequences(genome)
                self.clear_sequences(keep_proteins=True)

        self.calculate_gc_content()

        if custom_path == "":
            export_folder = self.path + "out_stats/"
        else:
            export_folder = custom_path
        if export_folder[-1] != "/":
            export_folder += "/"

        if export:
            system(f"mkdir -p {export_folder}")

        to_tally = ["coding_genes", "noncoding_genes", "CDSs_without_stop", "CDSs_with_stop"]

        self.stats = {"mean_transcripts" : [], "mean_exons" : [], "mean_exon_size" : [], "mean_gene_size" : [], "mean_intron_size" : [], "mean_CDS_size" : [], "mean_UTR_size" : [], "mean_transcript_size" : [], "mean_five_prime_UTR_size" : [], "mean_three_prime_UTR_size" : []}

        if genome != None:
            self.stats["mean_protein_size"] = []

        for key in to_tally:
            self.stats[key] = []

        for ft in self.features:
            self.stats[ft] = 0

        self.stats["five_prime_UTRs"] = 0
        self.stats["three_prime_UTRs"] = 0

        for ft in self.atypical_features:
            self.stats[ft.feature] = 0

        for ft in self.atypical_features:
            self.stats[ft.feature] += 1

        for genes in self.chrs.values():
            for g in genes.values():
                self.stats[g.feature] += 1
                self.stats["mean_transcripts"].append(len(g.transcripts))
                self.stats["mean_gene_size"].append(g.size)
                for t in g.transcripts.values():
                    self.stats[t.feature] += 1
                    if t.main:
                        if hasattr(t, "introns"):
                            for i in t.introns:
                                self.stats["mean_intron_size"].append(i.size)

                        if t.coding:
                            for c in t.CDSs.values():
                                self.stats[c.feature] += 1
                                if c.main:
                                    self.stats["mean_CDS_size"].append(c.size)
                                    if hasattr(c, "UTRs"):
                                        utr_5 = False
                                        utr_3 = False
                                        u_size = 0
                                        u5_size = 0
                                        u3_size = 0
                                        for u in c.UTRs:
                                            if u.prime == "3'":
                                                utr_3 = True
                                                u3_size += u.size
                                            else:
                                                utr_5 = True
                                                u5_size += u.size
                                            u_size += u.size
                                        if utr_5:
                                            self.stats["five_prime_UTRs"] += 1
                                        if utr_3:
                                            self.stats["three_prime_UTRs"] += 1
                                        self.stats["mean_UTR_size"].append(u_size)
                                        self.stats["mean_five_prime_UTR_size"].append(u5_size)
                                        self.stats["mean_three_prime_UTR_size"].append(u3_size)
                                    if genome != None:
                                        self.stats["mean_protein_size"].append(c.protein.size)
                            self.stats["coding_genes"].append(g.id)
                        else:
                            self.stats["noncoding_genes"].append(g.id)

                        for e in t.exons:
                            self.stats[e.feature] += 1
                            self.stats["mean_exon_size"].append(e.size)

                        self.stats["mean_exons"].append(len(t.exons))
                        
                        self.stats["mean_transcript_size"].append(t.size)

        self.stats["UTRs"] = self.stats["five_prime_UTRs"] + self.stats["three_prime_UTRs"]

        # anything with mean will be also plotted as distribution plots:
        if export:
            for key in self.stats:
                if "mean" in key:
                    tag = key.split("mean_")[1]
                    if tag[-1] != "s":
                        tag += "s"
                    barplot(self.stats[key], export_folder, f"{self.id}{self.feature_suffix}_{tag}", f"Distribution of {self.id} {tag}", max_x)   

        if genome != None:
            self.stats["CDSs_without_stop"] = []
            self.stats["CDSs_with_stop"] = []
            self.stats["intron_composition"] = set()
            intron_stats = {}
            for b in Intron.canonical_seqs:
                intron_stats[f"intron-exon boundary: {b}"] = 0
            intron_stats["other_intron_seqs"] = 0
            for genes in self.chrs.values():
                for g in genes.values():                    
                    for t in g.transcripts.values():
                        if t.main:
                            if t.coding:
                                for c in t.CDSs.values():
                                    if c.main:
                                        if c.protein != None:
                                            if not c.protein.end_stop:
                                                self.stats["CDSs_without_stop"].append(g.id)
                                            else:
                                                self.stats["CDSs_with_stop"].append(g.id)
                                            
                            for i in t.introns:
                                self.stats["intron_composition"].add(i.boundary)
                                if i.canonical:
                                    intron_stats[f"intron-exon boundary: {i.boundary}"] += 1
                                else:
                                    intron_stats["other_intron_seqs"] += 1
            self.stats["intron_composition"] = list(self.stats["intron_composition"])

            if export:
                labels = list(intron_stats.keys())
                mod_labels = []
                for l in labels:
                    if "other" in l:
                        mod_labels.append("other")
                    else:
                        mod_labels.append(l.split(": ")[1])

                values = list(intron_stats.values())
                pie_chart(mod_labels, values, export_folder, f"{self.id}{self.feature_suffix}_intron_composition", f"Intron composition of {self.id} annotation")
            self.stats.update(intron_stats)

        for key in self.stats:
            if "mean" in key:
                if len(self.stats[key]) > 0:
                    self.stats[key] = mean(self.stats[key])
                else:
                    self.stats[key] = 0
            #tallying
            elif key in to_tally:
                self.stats[key] = len(self.stats[key])

        if export:
            if genome != None:
                out_file = f"{self.id}{self.feature_suffix}_full_stats.csv"
            else:
                out_file = f"{self.id}{self.feature_suffix}_basic_stats.csv"

            warnings_df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in self.warnings.items()])).fillna('')
            error_df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in self.errors.items()])).fillna('')

            warnings_df.to_csv(f"{export_folder}{self.id}{self.feature_suffix}_warnings.csv", sep="\t", index=False)
            error_df.to_csv(f"{export_folder}{self.id}{self.feature_suffix}_errors.csv", sep="\t", index=False)

            f_out = open(f"{export_folder}{out_file}", "w", encoding="utf-8")
            f_out.write("")
            f_out.close()

            f_out = open(f"{export_folder}{out_file}", "a")
            x = -1
            for key, value in self.stats.items():
                x += 1
                value_temp = value
                if isinstance(value, list):
                    if len(value) > 200:
                        value_temp = value[:199]
                        value_temp.append("...")
                if x == 0:
                    f_out.write(f"{key}\t{value_temp}")
                else:
                    f_out.write(f"\n{key}\t{value_temp}")
            f_out.close()

        #print(f"\nUpdated stats for {self.id}")

    def detect_gene_overlaps(self, other:object=None, sort_processes:int=2, clear=True):
        """
        Detecting gene overlaps within the same annotation object or between
        annotation objects, provided they refer to the same genome.
        """
        start = time.time()
        # Check if stdout or stderr are redirected to files
        stdout_redirected = not sys.stdout.isatty()
        stderr_redirected = not sys.stderr.isatty()

        # Disable tqdm if stdout or stderr are redirected
        if stdout_redirected or stderr_redirected:
            disable = True
        else:
            disable = False

        if not self.sorted:
            self.sort_genes(processes=sort_processes)
        
        if other != None:
            if not other.sorted:
                other.sort_genes(processes=sort_processes)

        if clear:
            self.clear_overlaps()
            if other != None:
                other.clear_overlaps()

        if other != None:

            if self.genome == other.genome:
                
                if other.name in self.overlapped_annotations or self.name in other.overlapped_annotations:
                    print(f"Overlaps between {self.id} and {other.id} "
                           "annotations have already been detected, please "
                           "run 'self.clear_overlaps()' if you want to "
                           "recalculate them")
                else:
                    start = time.time()
                    progress_bar = tqdm(total=len(self.all_gene_ids.keys()), disable=disable,
                                bar_format=(
                    f'\033[1;91mWorking out overlaps between {self.id} and {other.id} annotations:\033[0m '
                    '{percentage:3.0f}%|'
                    f'\033[1;91m{{bar}}\033[0m| '
                    '{n}/{total} [{elapsed}<{remaining}]'))

                    if other.name not in self.overlapped_annotations:
                        self.overlapped_annotations.append(other.name)
                    if self.name not in other.overlapped_annotations:
                        other.overlapped_annotations.append(self.name)
                    for chr, genes in self.chrs.items():
                        for g1 in genes.values():
                            progress_bar.update(1)
                            found_overlap = False
                            if chr in other.chrs:
                                for g2 in other.chrs[chr].values():
                                    overlapping, overlap_bp = overlap(g1, g2)
                                    if overlapping:
                                        found_overlap = True
                                    elif found_overlap and g1.end < g2.start:
                                        break
                                    else:
                                        continue
                                    if g1.strand == g2.strand:
                                        gene_orientation = True
                                    else:
                                        gene_orientation = False

                                    gene_query_percent = (overlap_bp/(g1.size))*100                  
                                    gene_target_percent = (overlap_bp/(g2.size))*100

                                    target_exons = False
                                    query_exons = False
                                    best_exon_overlap = 0
                                    exon_orientation = False
                                    overlapping = False
                                    for t1 in g1.transcripts.values():
                                        if t1.exons != []:
                                            query_exons = True
                                        for t2 in g2.transcripts.values():
                                            if t2.exons != []:
                                                target_exons = True
                                            if t1.strand == t2.strand:
                                                exon_orientation = True
                                            overlap_exon_temp = 0
                                            for e1 in t1.exons:
                                                for e2 in t2.exons:
                                                    overlap_temp, overlap_bp = (
                                                                        overlap(e1, e2))
                                                    if overlap_temp:
                                                        overlap_exon_temp += overlap_bp
                                                        overlapping = True
                                            if overlap_exon_temp > best_exon_overlap:
                                                best_exon_overlap = overlap_exon_temp
                                                exon_query_size = t1.size
                                                exon_target_size = t2.size
                                    
                                    if target_exons and query_exons:
                                        exons_in_both = True
                                        if gene_orientation != exon_orientation:
                                            print(f"Warning: {self.id} query and {other.id} target have discrepancies in the orientation of gene and exons. Genes: {g1.id} and {g2.id}.")
                                        if overlapping:
                                            exon_query_percent = (best_exon_overlap/exon_query_size)*100
                                            exon_target_percent = (best_exon_overlap/exon_target_size)*100
                                        else:
                                            exon_query_percent = 0
                                            exon_target_percent = 0
                                    else:
                                        exons_in_both = False
                                        exon_orientation = None
                                        exon_query_percent = None
                                        exon_target_percent = None


                                    target_CDS = False
                                    query_CDS = False
                                    best_CDS_overlap = 0
                                    CDS_orientation = False
                                    overlapping = False
                                    for t1 in g1.transcripts.values():
                                        if t1.CDSs != {}:
                                            query_CDS = True
                                        for t2 in g2.transcripts.values():
                                            if t2.CDSs != {}:
                                                target_CDS = True
                                            if t1.strand == t2.strand:
                                                CDS_orientation = True
                                            for CDS1 in t1.CDSs.values():
                                                for CDS2 in t2.CDSs.values():
                                                    overlap_CDS_temp = 0
                                                    for c1 in CDS1.CDS_segments:
                                                        for c2 in CDS2.CDS_segments:
                                                            (overlap_temp,
                                                                overlap_bp) = (
                                                                        overlap(c1, c2))
                                                            if not overlap_temp:
                                                                continue
                                                            overlap_CDS_temp += (
                                                                            overlap_bp)
                                                            overlapping = True
                                                    if (overlap_CDS_temp > 
                                                        best_CDS_overlap):
                                                        best_CDS_overlap = (
                                                                    overlap_CDS_temp)
                                                        CDS_query_size = CDS1.size
                                                        CDS_target_size = CDS2.size
                                                
                                    if target_CDS and query_CDS:
                                        CDSs_in_both = True
                                        if gene_orientation != CDS_orientation:
                                            print(f"Warning: {self.id} query and {other.id} target have discrepancies in the orientation of gene and CDS. Genes: {g1.id} and {g2.id}.")
                                        if overlapping:
                                            CDS_query_percent = ((best_CDS_overlap
                                                                /CDS_query_size)*100)
                                            CDS_target_percent = ((best_CDS_overlap
                                                                /CDS_target_size)*100)
                                        else:
                                            CDS_query_percent = 0
                                            CDS_target_percent = 0
                                    else:
                                        CDSs_in_both = False
                                        CDS_orientation = None
                                        CDS_query_percent = None
                                        CDS_target_percent = None


                                    protein_query_percent = None
                                    protein_target_percent = None
                                    if CDS_query_percent != None and CDS_query_percent != 0:
                                        if CDS_orientation:
                                            target_protein = False
                                            query_protein = False
                                            best_protein_overlap = 0
                                            overlapping = False
                                            for t1 in g1.transcripts.values():
                                                if t1.CDSs != {}:
                                                    query_protein = True
                                                for t2 in g2.transcripts.values():
                                                    if t2.CDSs != {}:
                                                        target_protein = True
                                                    for CDS1 in t1.CDSs.values():
                                                        for CDS2 in t2.CDSs.values():
                                                            overlap_protein_temp = 0
                                                            for c1 in CDS1.CDS_segments:
                                                                for c2 in CDS2.CDS_segments:
                                                                    (overlap_temp,
                                                                        overlap_bp) = overlap(c1,
                                                                                                c2)
                                                                    if not overlap_temp:
                                                                        continue
                                                                    if c1.frame != c2.frame:
                                                                        continue
                                                                    
                                                                    overlap_protein_temp += overlap_bp
                                                                    overlapping = True
                                                            if overlap_protein_temp > best_protein_overlap:
                                                                best_protein_overlap = overlap_protein_temp
                                                                protein_query_size = CDS1.size
                                                                protein_target_size = CDS2.size
                                                        
                                            if target_protein and query_protein:
                                                if overlapping:
                                                    protein_query_percent = ((best_protein_overlap
                                                                            /protein_query_size)*100)
                                                    protein_target_percent = ((best_protein_overlap
                                                                            /protein_target_size)*100)
                                                else:
                                                    protein_query_percent = 0
                                                    protein_target_percent = 0

                                    g1.overlaps["other"].append(OverlapHit(g2.id, 
                                                                    other.name,
                                                                    gene_orientation,
                                                                    gene_query_percent,
                                                                    gene_target_percent,
                                                                    exons_in_both,
                                                                    exon_query_percent,
                                                                    exon_target_percent,
                                                                    CDSs_in_both,
                                                                    CDS_query_percent,
                                                                    CDS_target_percent,
                                                                    protein_query_percent,
                                                                    protein_target_percent,
                                                                g2.conserved_synteny))
                                    g2.overlaps["other"].append(OverlapHit(g1.id,
                                                                    self.name,
                                                                    gene_orientation,
                                                                    gene_target_percent,
                                                                    gene_query_percent,
                                                                    exons_in_both,
                                                                    exon_target_percent,
                                                                    exon_query_percent,
                                                                    CDSs_in_both,
                                                                    CDS_target_percent,
                                                                    CDS_query_percent,
                                                                    protein_query_percent,
                                                                    protein_target_percent,
                                                                g1.conserved_synteny))
                    self.add_aliases()
                    other.add_aliases()
                    now = time.time()
                    lapse = now - start
                    progress_bar.close()
                    print(f"Detecting overlaps between {other.id} and {self.id} annotations took {round(lapse/60, 1)} minutes")
            else:
                print(f"Did not generate overlaps between {other.id} and {self.id} annotations as they are associated to different genomes")

        else:
            if self.self_overlapping != []:
                print("There are already detected 'self' gene overlaps, please run 'self.clear_overlaps()' if you want to recalculate them")
            else:
                progress_bar = tqdm(total=len(self.all_gene_ids.keys()), disable=disable,
                            bar_format=(
                f'\033[1;91mWorking out overlaps within {self.id} annotation:\033[0m '
                '{percentage:3.0f}%|'
                f'\033[1;91m{{bar}}\033[0m| '
                '{n}/{total} [{elapsed}<{remaining}]'))

                # making sure self overlaps are not added twice
                start = time.time()
                self.self_overlapping = set(self.self_overlapping)
                for chr, genes in self.chrs.items():
                    gl = list(genes.keys())[1:]
                    for g1 in genes.values():
                        progress_bar.update(1)
                        found_overlap = False
                        for gl_id in gl:
                            g2 = genes[gl_id]
                            if g1.id == g2.id:
                                continue
                            overlapping, overlap_bp = overlap(g1, g2)

                            if overlapping:
                                self.self_overlapping.add(g1.id)
                                self.self_overlapping.add(g2.id)
                                found_overlap = True

                            elif found_overlap and g1.end < g2.start:
                                break
                            else:
                                continue

                            if g1.strand == g2.strand:
                                gene_orientation = True
                            else:
                                gene_orientation = False

                            gene_query_percent = (overlap_bp/(g1.size))*100                  
                            gene_target_percent = (overlap_bp/(g2.size))*100

                            target_exons = False
                            query_exons = False
                            best_exon_overlap = 0
                            exon_orientation = False
                            overlapping = False
                            for t1 in g1.transcripts.values():
                                if t1.exons != []:
                                    query_exons = True
                                for t2 in g2.transcripts.values():
                                    if t2.exons != []:
                                        target_exons = True
                                    if t1.strand == t2.strand:
                                        exon_orientation = True
                                    overlap_exon_temp = 0
                                    for e1 in t1.exons:
                                        for e2 in t2.exons:
                                            overlap_temp, overlap_bp = overlap(e1,
                                                                            e2)
                                            if overlap_temp:
                                                overlap_exon_temp += overlap_bp
                                                overlapping = True
                                    if overlap_exon_temp > best_exon_overlap:
                                        best_exon_overlap = overlap_exon_temp
                                        exon_query_size = t1.size
                                        exon_target_size = t2.size
                            
                            if target_exons and query_exons:
                                exons_in_both = True
                                if gene_orientation != exon_orientation:
                                    print(f"Warning: {self.id} query and target have discrepancies in the orientation of gene and exons. Genes: {g1.id} and {g2.id}")
                                if overlapping:
                                    exon_query_percent = ((best_exon_overlap
                                                            /exon_query_size)*100)
                                    exon_target_percent = ((best_exon_overlap
                                                            /exon_target_size)*100)
                                else:
                                    exon_query_percent = 0
                                    exon_target_percent = 0
                            else:
                                exons_in_both = False
                                exon_orientation = None
                                exon_query_percent = None
                                exon_target_percent = None

                            target_CDS = False
                            query_CDS = False
                            best_CDS_overlap = 0
                            CDS_orientation = False
                            overlapping = False
                            for t1 in g1.transcripts.values():
                                if t1.CDSs != {}:
                                    query_CDS = True
                                for t2 in g2.transcripts.values():
                                    if t2.CDSs != {}:
                                        target_CDS = True
                                    if t1.strand == t2.strand:
                                        CDS_orientation = True
                                    for CDS1 in t1.CDSs.values():
                                        for CDS2 in t2.CDSs.values():
                                            overlap_CDS_temp = 0
                                            for c1 in CDS1.CDS_segments:
                                                for c2 in CDS2.CDS_segments:
                                                    (overlap_temp,
                                                        overlap_bp) = overlap(c1,
                                                                                c2)
                                                    if not overlap_temp:
                                                        continue
                                                    overlap_CDS_temp += overlap_bp
                                                    overlapping = True
                                            if overlap_CDS_temp > best_CDS_overlap:
                                                best_CDS_overlap = overlap_CDS_temp
                                                CDS_query_size = CDS1.size
                                                CDS_target_size = CDS2.size
                                        
                            if target_CDS and query_CDS:
                                CDSs_in_both = True
                                if gene_orientation != CDS_orientation:
                                    print(f"Error: {self.id} query and target have discrepancies"
                                            " in the orientation of gene and CDS."
                                            f"genes: {g1.id} and {g2.id}" 
                                            "DO NOT CONTINUE! -> fix the problem!")
                                if overlapping:
                                    CDS_query_percent = ((best_CDS_overlap
                                                            /CDS_query_size)*100)
                                    CDS_target_percent = ((best_CDS_overlap
                                                            /CDS_target_size)*100)
                                else:
                                    CDS_query_percent = 0
                                    CDS_target_percent = 0
                            else:
                                CDSs_in_both = False
                                CDS_orientation = None
                                CDS_query_percent = None
                                CDS_target_percent = None

                            if CDSs_in_both:
                                protein_query_percent = 0
                                protein_target_percent = 0
                                if CDS_query_percent != None and CDS_query_percent != 0:
                                    if CDS_orientation:
                                        target_protein = False
                                        query_protein = False
                                        best_protein_overlap = 0
                                        overlapping = False
                                        for t1 in g1.transcripts.values():
                                            if t1.CDSs != {}:
                                                query_protein = True
                                            for t2 in g2.transcripts.values():
                                                if t2.CDSs != {}:
                                                    target_protein = True
                                                for CDS1 in t1.CDSs.values():
                                                    for CDS2 in t2.CDSs.values():
                                                        overlap_protein_temp = 0
                                                        for c1 in CDS1.CDS_segments:
                                                            for c2 in CDS2.CDS_segments:
                                                                (overlap_temp,
                                                                    overlap_bp) = overlap(c1,
                                                                                            c2)
                                                                if not overlap_temp:
                                                                    continue
                                                                if c1.frame != c2.frame:
                                                                    continue
                                                                
                                                                overlap_protein_temp += overlap_bp
                                                                overlapping = True
                                                        if overlap_protein_temp > best_protein_overlap:
                                                            best_protein_overlap = overlap_protein_temp
                                                            protein_query_size = CDS1.size
                                                            protein_target_size = CDS2.size
                                                    
                                        if target_protein and query_protein:
                                            if overlapping:
                                                protein_query_percent = ((best_protein_overlap
                                                                        /protein_query_size)*100)
                                                protein_target_percent = ((best_protein_overlap
                                                                        /protein_target_size)*100)
                            else:

                                protein_query_percent = None
                                protein_target_percent = None
                            
                            
                            g1.overlaps["self"].append(OverlapHit(g2.id, self.name,
                                                                gene_orientation,
                                                                gene_query_percent,
                                                                gene_target_percent,
                                                                exons_in_both,
                                                                exon_query_percent,
                                                                exon_target_percent,
                                                                CDSs_in_both,
                                                                CDS_query_percent,
                                                                CDS_target_percent,
                                                                protein_query_percent,
                                                                protein_target_percent,
                                                            g2.conserved_synteny))
                            g2.overlaps["self"].append(OverlapHit(g1.id, self.name,
                                                                gene_orientation,
                                                                gene_target_percent,
                                                                gene_query_percent,
                                                                exons_in_both,
                                                                exon_target_percent,
                                                                exon_query_percent,
                                                                CDSs_in_both,
                                                                CDS_target_percent,
                                                                CDS_query_percent,
                                                                protein_query_percent,
                                                                protein_target_percent,
                                                            g1.conserved_synteny))
                        try:
                            gl.remove(g1.id)
                        except:
                            pass
                self.self_overlapping = list(self.self_overlapping)
                progress_bar.close()
                now = time.time()
                lapse = now - start
                print(f"\nDetecting gene overlaps within the {self.id} "
                      f"annotation took {round(lapse/60, 1)} minutes\n")
                print(f"\nThere are {len(self.self_overlapping)} genes "
                      f"overlapping with other genes in {self.id} annotation\n")
                self.add_qualitative_info_to_overlaps()

    def overlaps_as_networks(self, self_mode:bool=True):
        self.overlap_networks = {}
        for chr, genes in self.chrs.items():
            G = nx.Graph()
            for g in genes.values():
                if self_mode:
                    overlaps = g.overlaps["self"]
                else:
                    overlaps = g.overlaps["other"]
                for o in overlaps:
                    G.add_edge(g.id, o.id)
            self.overlap_networks[chr] = list(nx.connected_components(G))

    def alternative_remove_redundancy(self):
        nodes = self.overlap_networks[chr][0].nodes()
        print("Nodes in the graph:")
        for node in nodes:
            print(node)

        # Find articulation points (connector nodes)
        connector_nodes = list(nx.articulation_points(self.overlap_networks[chr][0]))

        # Remove connector nodes from the graph
        for node in connector_nodes:
            self.overlap_networks[chr][0].remove_node(node)

    def add_aliases(self, stringent:bool=True):
        if stringent:
            acceptable_values = [6, 7, 8, 9, 10, 11]
        else:
            acceptable_values = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
        for genes in self.chrs.values():
            for g in genes.values():
                for name, hits in g.overlaps.items():
                    if name == "other":
                        for hit in hits:
                            if hit.score in acceptable_values:
                                if hit.id not in g.aliases:
                                    g.aliases.append(hit.id)

    def export_equivalences(self, custom_path:str="", stringent:bool=True, verbose:bool=True, return_df:bool=True, NAs:bool=True, export_csv:bool=False, export_self:bool=False):
        start = time.time()
        if export_self:
            export = "self"
            export_tag = "self_"
        else:
            export = "other"
            export_tag = ""
        tag = f"{export_tag}{self.id}{self.feature_suffix}_equivalences"
        if stringent:
            acceptable_values = [6, 7, 8, 9, 10, 11]
            tag += "_filtered"
        else:
            acceptable_values = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]        
        
        if custom_path == "":
            export_folder = self.path + "overlaps/"
        else:
            export_folder = custom_path
            
        if export_folder[-1] != "/":
            export_folder += "/"
        system(f"mkdir -p {export_folder}")

        equivalences = {"query_id": [], "target_id": []}
        if verbose:
            equivalences["query_synteny_conserved"] = []
            equivalences["target_synteny_conserved"] = []
            equivalences["same_strand"] = []
            equivalences["min_gene_percent"] = []
            equivalences["min_exon_percent"] = []
            equivalences["min_CDS_percent"] = []
        if not verbose:
            tag += "_simple"
        equivalences["query_origin"] = []            
        equivalences["target_origin"] = []
        equivalences["overlap_score"] = []
        
        eq_df = pd.DataFrame(data=equivalences)
        for genes in self.chrs.values():
            for g in genes.values():
                for name, hits in g.overlaps.items():
                    if name == export:
                        for hit in hits:
                            if hit.score not in acceptable_values:
                                continue
                            index = len(eq_df.index)
                            if eq_df.empty:
                                eq_df.loc[0] = pd.NA
                            eq_df.loc[index, "query_id"] = g.id
                            eq_df.loc[index, "target_id"] = hit.id
    
                            
                            if verbose:      
                                eq_df.loc[index, "same_strand"] = hit.orientation
                                if self.liftoff:
                                    eq_df.loc[index, "query_synteny_conserved"] = g.conserved_synteny
                                else:
                                    eq_df.loc[index, "query_synteny_conserved"] = pd.NA
                                eq_df.loc[index, "target_synteny_conserved"] = hit.target_synteny_conserved

                                eq_df.loc[index, "min_gene_percent"] = hit.min_gene_percent
                                eq_df.loc[index, "min_exon_percent"] = hit.min_exon_percent
                                eq_df.loc[index, "min_CDS_percent"] = hit.min_CDS_percent
                            eq_df.loc[index, "query_origin"] = self.name
                            eq_df.loc[index, "target_origin"] = hit.origin
                            eq_df.loc[index, "overlap_score"] = hit.score

        if NAs:
            tag += "_justQueryNAs"
            overlapping_genes = eq_df["query_id"].to_list()
            # adding target genes with no equivalence:
            for genes in self.chrs.values():
                for g in genes.values():
                    if g.id not in overlapping_genes:
                        index = len(eq_df.index)
                        eq_df.loc[index] = pd.NA
                        eq_df.loc[index, "query_id"] = g.id
                        eq_df.loc[index, "query_origin"] = self.name
                        eq_df.loc[index, "overlap_score"] = 0
            # Only for liftoff annotations
            for g_id in self.unmapped:
                if g_id not in overlapping_genes:
                    index = len(eq_df.index)
                    eq_df.loc[index] = pd.NA
                    eq_df.loc[index, "query_id"] = g_id
                    eq_df.loc[index, "query_origin"] = self.name
        else:
            tag += "_noNAs"

        # There should not be any duplicates but removing them just in case
        eq_df.drop_duplicates(inplace=True)
        column_sort_order = ["query_origin", "target_origin", "overlap_score", "query_synteny_conserved", "target_synteny_conserved", "query_id", "target_id"]
        ascending = [True, True, False, False, False, True, True]
        eq_df.sort_values(by=column_sort_order, ascending=ascending, inplace=True)
        eq_df.reset_index(drop=True, inplace=True)
        
        if export_csv:
            eq_df.to_csv(f"{export_folder}{tag}.csv", sep="\t", index=False, na_rep="NA")

        now = time.time()
        lapse = now - start
        print(f"\nExporting {self.id} equivalences to the following queries {self.overlapped_annotations} took {round(lapse/60, 1)} minutes")
        
        if return_df:
            return eq_df
                  
    def clear_overlaps(self, keep_self=False, keep_other=False):
        if not keep_self:
            self.self_overlapping = []
            for genes in self.chrs.values():
                for g in genes.values():
                    g.overlaps["self"] = []
        if not keep_other:
            self.overlapped_annotations = []
            for genes in self.chrs.values():
                for g in genes.values():
                    g.overlaps["other"] = []
    
    def clear_aliases(self):
        for genes in self.chrs.values():
            for g in genes.values():
                g.aliases = []

    def export_lengths(self, just_genes:bool=True, custom_path:str=""):
        """
        Exports gene lengths as total exon length of main transcript.
        """
        out = ""
        if custom_path == "":
            export_folder = self.path + "gene_lengths/"
        else:
            export_folder = custom_path
        if export_folder[-1] != "/":
            export_folder += "/"
        system(f"mkdir -p {export_folder}")
        if not just_genes:
            for genes in self.chrs.values():
                for g in genes.values():
                    for t in g.transcripts.values():
                        if t.main:
                            out += f"{g.id}\t{t.size}\n"
            tag = f"{self.id}{self.feature_suffix}_gene_lengths_all.tsv"
        else:
            for genes in self.chrs.values():
                for g in genes.values():
                    if not g.pseudogene and not g.transposable:
                        for t in g.transcripts.values():
                            if t.main:
                                out += f"{g.id}\t{t.size}\n"
            tag = f"{self.id}{self.feature_suffix}_gene_lengths.tsv"
        f_out = open(f"{export_folder}{tag}", "w", encoding="utf-8")
        f_out.write(out)
        f_out.close()

    def export_coordinates(self, custom_path:str="", lengths:bool=False):
        """
        Exports gene coordinates.
        """
        if lengths:
            out = "gene_id\tchromosome\tgene_start\tgene_end\tgene_length\n"
        else:
            out = "gene_id\tchromosome\tgene_start\tgene_end\n"
        if custom_path == "":
            export_folder = self.path + "gene_coordinates/"
        else:
            export_folder = custom_path
        if export_folder[-1] != "/":
            export_folder += "/"
        system(f"mkdir -p {export_folder}")
        if lengths:
            for chrom, genes in self.chrs.items():
                for g in genes.values():
                    out += f"{g.id}\t{chrom}\t{g.start}\t{g.end}\t{g.size}\n"    
        else:
            for chrom, genes in self.chrs.items():
                for g in genes.values():
                    out += f"{g.id}\t{chrom}\t{g.start}\t{g.end}\n"
        if self.dapmod:
            tag = f"{self.id}{self.feature_suffix}_gene_coordinates_dapmod.tsv"
        else:
            tag = f"{self.id}{self.feature_suffix}_gene_coordinates.tsv"
        f_out = open(f"{export_folder}{tag}", "w", encoding="utf-8")
        f_out.write(out)
        f_out.close()

    def CDS_to_CDS_segment_ids(self, extra_attributes:bool=False):
        for genes in self.chrs.values():
            for g in genes.values():
                for t in g.transcripts.values():
                    for c in t.CDSs.values():
                        for x, cs in enumerate(c.CDS_segments):
                            cs.id = f"{cs.id}_{x+1}"

        self.update_attributes(extra_attributes=extra_attributes)

    def export_gff(self, custom_path:str="", tag:str=".gff3", skip_atypical_fts:bool=True, main_only:bool=False, UTRs:bool=False, just_genes:bool=False, no_1bp_features:bool=False):

        # Check if stdout or stderr are redirected to files
        stdout_redirected = not sys.stdout.isatty()
        stderr_redirected = not sys.stderr.isatty()

        # Disable tqdm if stdout or stderr are redirected
        if stdout_redirected or stderr_redirected:
            disable = True
        else:
            disable = False

        progress_bar = tqdm(total=len(self.all_gene_ids.keys()), disable=disable,
                                bar_format=(
                    f'\033[38;2;46;204;113m\033[1mExporting gff for {self.id} annotation:\033[0m '
                    '{percentage:3.0f}%|'
                    f'\033[38;2;46;204;113m{{bar}}\033[0m| '
                    '{n}/{total} [{elapsed}<{remaining}]'))

        out = ""
        if custom_path == "":
            export_folder = self.path + "out_gffs/"
        else:
            export_folder = custom_path
        if export_folder[-1] != "/":
            export_folder += "/"

        system(f"mkdir -p {export_folder}")

        if self.clean or self.dapmod:
            out += "##gff-version 3\n"
        else:
            out += self.gff_header

        for genes in self.chrs.values():
            progress_bar.update(len(genes))
            for g in genes.values():
                
                if just_genes:
                    continue
                
                if no_1bp_features:
                    gene_1bp_feature = False
                    for t in g.transcripts.values():
                        for e in t.exons:
                            if e.size == 1:
                                gene_1bp_feature = True
                        for c in t.CDSs.values():
                            for cs in c.CDS_segments:
                                if cs.size == 1:
                                    gene_1bp_feature = True

                    if gene_1bp_feature:
                        continue

                out += g.print_gff()

                temp_subfeatures = []
                for t in g.transcripts.values():
                    if main_only:
                        if not t.main:
                            continue
                    out += t.print_gff()
                    for e in t.exons:
                        add = True
                        for ts in temp_subfeatures:
                            if e.almost_equal(ts):
                                add = False
                        if add:
                            temp_subfeatures.append(e)
                    for c in t.CDSs.values():
                        if main_only:
                            if not c.main:
                                continue
                        for c_seg in c.CDS_segments:
                            add = True
                            for ts in temp_subfeatures:
                                if c_seg == ts:
                                    add = False
                            if add:
                                temp_subfeatures.append(c_seg)
                        if UTRs:
                            if hasattr(c, "UTRs"):
                                for u in c.UTRs:
                                    add = True
                                    for ts in temp_subfeatures:
                                        if u == ts:
                                            add = False
                                    if add:
                                        temp_subfeatures.append(u)
                temp_subfeatures.sort()
                for ts in temp_subfeatures:
                    out += ts.print_gff()
                out += "###\n"

        progress_bar.close()

        if not skip_atypical_fts:
            if not just_genes:
                if self.atypical_features != []:
                    print("Export of atypical features is not fully supported yet "
                        "these will just be added at the end of the gff regardless "
                        "of coordinate or chromosome/scaffold order")
                    for ft in self.atypical_features:
                        out += ft.print_gff()

        output_suffix = ""

        if just_genes:
            output_suffix += "_just_genes"
        if no_1bp_features:
            output_suffix += "_for_lifton"

        if tag == ".gff3":
            tag = f"{self.id}{self.suffix}{output_suffix}{tag}"


        print(f"Exporting {self.id} gff with tag='{tag}' which is dapfit={self.dapfit} and dapmod={self.dapmod} and combined={self.combined}.")         

        f_out = open(f"{export_folder}{tag}", "w", encoding="utf-8")
        f_out.write(out)
        f_out.close()
        
    def merge(self, other:object, ignore_overlaps:bool=True, exon_overlap_threshold:float=100, gene_overlap_threshold:float=100):
        """
        Priority is given to self annotation
        """
        start = time.time()
        self.update()
        other.update()

        if not ignore_overlaps:
            self.detect_gene_overlaps(other)

        # Check if stdout or stderr are redirected to files
        stdout_redirected = not sys.stdout.isatty()
        stderr_redirected = not sys.stderr.isatty()

        # Disable tqdm if stdout or stderr are redirected
        if stdout_redirected or stderr_redirected:
            disable = True
        else:
            disable = False

        progress_bar = tqdm(total=len(other.all_gene_ids.keys()), disable=disable,
                                bar_format=(
                    f'\033[38;2;46;204;113m\033[1mMerging {other.id} and {self.id} annotations:\033[0m '
                    '{percentage:3.0f}%|'
                    f'\033[38;2;46;204;113m{{bar}}\033[0m| '
                    '{n}/{total} [{elapsed}<{remaining}]'))
        if self.merged:
            if "_..._" in self.name:
                first_name = self.name.split("_..._")[0]
            else:
                first_name = self.name.split("_merged_")[0]
            self.name = f"{first_name}_..._{other.name}"
        else:
            self.name = f"{self.name}_merged_{other.name}"
        self.id = f"{self.name}_on_{self.genome}"
        count = 0
        if ignore_overlaps:
            for chr, genes in other.chrs.items():
                if chr not in self.chrs:
                    self.chrs[chr] = {}
                for g in genes.values():
                    progress_bar.update(1)
                    temp_id = g.id
                    while temp_id in self.all_gene_ids:
                        count += 1
                        temp_id = f"{g.id}_{count}"
                    self.chrs[chr][temp_id] = g.copy()
                    self.chrs[chr][temp_id].id = temp_id

        else:
            for chr, genes in other.chrs.items():
                if chr not in self.chrs:
                    self.chrs[chr] = {}
                for g in genes.values():
                    progress_bar.update(1)
                    count = 0
                    if g.overlaps["other"] == []:
                        temp_id = g.id
                        while temp_id in self.all_gene_ids:
                            count += 1
                            temp_id = f"{g.id}_{count}"
                        self.chrs[chr][temp_id] = g.copy()
                        self.chrs[chr][temp_id].id = temp_id
                    else:
                        exon_scores = [0]
                        gene_scores = [0]
                        for o in g.overlaps["other"]:
                            if o.min_exon_percent != None:
                                exon_scores.append(o.min_exon_percent)
                            if o.min_gene_percent != None:
                                gene_scores.append(o.min_gene_percent)
                        # some genes may not have annotated and hence exon_overlap_threshold cannot be taken into account
                        if o.exons_in_both:
                            if max(exon_scores) <= exon_overlap_threshold and max(gene_scores) <= gene_overlap_threshold:
                                temp_id = g.id
                                while temp_id in self.all_gene_ids:
                                    count += 1
                                    temp_id = f"{g.id}_{count}"
                                self.chrs[chr][temp_id] = g.copy()
                                self.chrs[chr][temp_id].id = temp_id
                        elif max(gene_scores) <= gene_overlap_threshold:
                            temp_id = g.id
                            while temp_id in self.all_gene_ids:
                                count += 1
                                temp_id = f"{g.id}_{count}"
                            self.chrs[chr][temp_id] = g.copy()
                            self.chrs[chr][temp_id].id = temp_id

        self.merged = True
        self.sorted = False
        self.generated_all_sequences = False
        self.contains_all_sequences = False        
        self.generated_CDS_sequences = False
        self.contains_CDS_sequences = False
        self.generated_protein_sequences = False
        self.contains_protein_sequences = False
        progress_bar.close()
        self.rename_ids(from_merge=True)
        self.update()
        now = time.time()
        lapse = now - start
        print(f"\nMerging {self.id} and {other.id} annotations took {round(lapse/60, 1)} minutes")

    def remove_wrongly_assigned_exons(self):
        for genes in self.chrs.values():
            for g in genes.values():
                for t in g.transcripts.values():
                    new_e = []
                    for e in t.exons:
                        if e.strand == t.strand:
                            new_e.append(e)
                    t.exons = new_e.copy()

        self.remove_transcripts_with_no_exons()
        self.remove_genes_with_no_transcripts()

    def remove_transcripts_with_no_exons(self):
        transcripts_to_remove = []
        for chrom, genes in self.chrs.items():
            for g in genes.values():
                for t in g.transcripts.values():
                    if t.exons == []:
                        transcripts_to_remove.append((chrom, g.id, t.id))
        for chrom, g_id, t_id in transcripts_to_remove:
            del self.chrs[chrom][g_id].transcripts[t_id]
            print(f"{t_id} Warning: transcript with no exons whichcould have been removed because of strand inconsistencies")
            self.warnings["transcript_with_no_exons"].append(t_id)

    def remove_genes_with_no_transcripts(self):
        genes_to_remove = []
        for chrom, genes in self.chrs.items():
            for g in genes.values():
                if g.transcripts == {}:
                    genes_to_remove.append((chrom, g.id))
        for chrom, g_id in genes_to_remove:
            del self.chrs[chrom][g_id]
            print(f"{g_id} Warning: gene with no transcripts which could have been removed because of strand inconsistencies of its exons")
            self.warnings["gene_with_no_transcripts"].append(g_id)

    def remove_missing_transcript_parent_references(self, extra_attributes=False):
        #print(f"Removing missing transcript parent references for {self.id} annotation.")
        self.remove_wrongly_assigned_exons()

        for genes in self.chrs.values():
            for g in genes.values():
                for t in g.transcripts.values():
                    for e in t.exons:
                        new_parents = []
                        for p in e.parents:
                            if p in self.all_transcript_ids:
                                if g.transcripts[p].strand == e.strand:
                                    new_parents.append(p)
                        e.parents = new_parents.copy()
                    for c in t.CDSs.values():
                        new_parents = []
                        for p in c.parents:
                            if p in self.all_transcript_ids:
                                new_parents.append(p)
                        c.parents = new_parents.copy()
                        for cs in c.CDS_segments:
                            new_parents = []
                            for p in cs.parents:
                                if p in self.all_transcript_ids:
                                    new_parents.append(p)
                            cs.parents = new_parents.copy()

        self.update_attributes(extra_attributes=extra_attributes)
        #print(f"Removed missing transcript parent references for {self.id} annotation.")

    def rework_CDSs(self, genome:object=None, low_memory:bool=True):
        start = time.time()
        if low_memory:
            self.clear_sequences()
        # Check if stdout or stderr are redirected to files
        stdout_redirected = not sys.stdout.isatty()
        stderr_redirected = not sys.stderr.isatty()

        # Disable tqdm if stdout or stderr are redirected
        if stdout_redirected or stderr_redirected:
            disable = True
        else:
            disable = False
        progress_bar = tqdm(total=len(self.all_gene_ids.keys()), disable=disable,
                                bar_format=(
                    f'\033[1;91mReworking {self.id} CDSs:\033[0m '
                    '{percentage:3.0f}%|'
                    f'\033[1;91m{{bar}}\033[0m| '
                    '{n}/{total} [{elapsed}<{remaining}]'))
    
        for genes in self.chrs.values():
            for g in genes.values():
                progress_bar.update(1)
                for t in g.transcripts.values():
                    t.generate_sequence(genome, low_memory)
                    t.generate_best_protein(genome)
                    t.generate_CDSs_based_on_ORF(low_memory)
                    for c in t.CDSs.values():
                        c.generate_sequence(genome, low_memory)
                    t.exon_update()
                    if t.coding_ratio < 0.80:
                        t.generate_sequence(genome, low_memory)
                        t.generate_best_protein(genome, must_have_stop=False)
                        t.generate_CDSs_based_on_ORF(low_memory)
                        for c in t.CDSs.values():
                            c.generate_sequence(genome, low_memory)
                    t.exon_update()

        progress_bar.close()
        self.update()
        now = time.time()
        lapse = now - start
        print(f"\nReworking CDSs for {self.id} took {round(lapse/60, 1)} "
               "minutes")

    def update_gene_and_transcript_list(self):
                # Check if stdout or stderr are redirected to files
        stdout_redirected = not sys.stdout.isatty()
        stderr_redirected = not sys.stderr.isatty()

        # Disable tqdm if stdout or stderr are redirected
        if stdout_redirected or stderr_redirected:
            disable = True
        else:
            disable = False

        total = 0
        for genes in self.chrs.values():
            total += len(genes.keys())

        progress_bar = tqdm(total=total, disable=disable,
                                bar_format=(
                    f'\033[1;95mUpdating {self.id} annotation gene and transcript lists:\033[0m '
                    '{percentage:3.0f}%|'
                    f'\033[1;95m{{bar}}\033[0m| '
                    '{n}/{total} [{elapsed}<{remaining}]'))
        self.all_gene_ids = {}
        self.all_transcript_ids = {}
        for chr, genes in self.chrs.items():
            for g in genes.values():
                progress_bar.update(1)
                self.all_gene_ids[g.id] = chr
                for t in g.transcripts.values():
                    self.all_transcript_ids[t.id] = chr
        progress_bar.close()

    def make_alternative_transcripts_into_genes(self):
        new_chrs = {}
        for chrom, genes in self.chrs.items():
            new_genes = {}
            for g in genes.values():
                new_transcripts = {}
                for t in g.transcripts.values():
                    if t.main:
                        new_transcripts[t.id] = t.copy()
                    else:
                        g_new = g.copy()
                        count = 0
                        while g_new.id in self.all_gene_ids:
                            count += 1
                            g_new.id = f"{g.id}_{count}"
                        
                        g_new.transcripts = {t.id : t.copy()}
                        g_new.start = t.start
                        g_new.end = t.end
                        new_genes[g_new.id] = g_new.copy()
                        self.all_gene_ids[g_new.id] = chrom

                g.transcripts = new_transcripts.copy()
            new_chrs[chrom] = new_genes.copy()

        for nc in new_chrs:
            for g in new_chrs[nc].values():
                self.chrs[nc][g.id] = g.copy()

        del new_chrs

        self.update_gene_and_transcript_list()
        self.update()

    def rename_ids(self, minimal:bool=False, prefix:str="", suffix:str="", spacer:int=100, from_merge:bool=False, extra_attributes:bool=False, remove_point_suffix:bool=False, correspondences:bool=False, custom_path:str="", cs_segment_ids:bool=False):
        """
        If a prefix is provided, everything is renamed. Otherwise just transcripts and subfeatures.
        Minimal mode just changes CDS or transcript IDs that do not contain part of the gene id.
        """

        if custom_path == "":
            export_folder = self.path + "out_gffs/"
        else:
            export_folder = custom_path
        if export_folder[-1] != "/":
            export_folder += "/"


        if from_merge:
            self.update_gene_and_transcript_list()

        for genes in self.chrs.values():
            for g in genes.values():
                for t in g.transcripts.values():
                    t.update()
                g.update()

        start = time.time()
        # Check if stdout or stderr are redirected to files
        stdout_redirected = not sys.stdout.isatty()
        stderr_redirected = not sys.stderr.isatty()

        # Disable tqdm if stdout or stderr are redirected
        if stdout_redirected or stderr_redirected:
            disable = True
        else:
            disable = False
        changed = set()
        progress_bar = tqdm(total=len(self.all_gene_ids.keys()), disable=disable,
                                bar_format=(
                    f'\033[38;2;156;42;42m\033[1mRenaming {self.id} Gene Ids:\033[0m '
                    '{percentage:3.0f}%|'
                    f'\033[38;2;156;42;42m{{bar}}\033[0m| '
                    '{n}/{total} [{elapsed}<{remaining}]'))
        if not minimal:
            if not prefix:
                for genes in self.chrs.values():
                    for g in genes.values():
                        progress_bar.update(1)
                        changed.add("transcripts")
                        changed.add("exons")
                        changed.add("CDSs")
                        changed.add("CDS_segments")

                        if g.id.endswith("_gene"):
                            g_prefix = g.id[:-5]
                            changed.add("genes")
                        elif g.id.endswith("gene"):
                            g_prefix = g.id[:-4]
                            changed.add("genes")
                        elif g.id.startswith("gene:"):
                            g_prefix = g.id[5:]
                            changed.add("genes")
                        elif g.id.startswith("gene-"):
                            g_prefix = g.id[5:]  
                        elif g.id.startswith("gene"):
                            g_prefix = g.id[4:]
                            changed.add("genes")
                        else:
                            g_prefix = g.id

                        if remove_point_suffix:
                            g.id = g.id.split(".")[0]
                            g_prefix = g.id

                        g.id = g_prefix
                        
                        t_count = 1
                        tmains = 0

                        for t in g.transcripts.values():
                            t.parents = [g_prefix]
                            if t.main:
                                t.id = f"{g_prefix}_t001"
                                tmains += 1
                            else:
                                t_count += 1
                                t_count_s = "{:03d}".format(t_count)
                                t.id = f"{g_prefix}_t{t_count_s}"
                            c_count = 1
                            # CDSs renamed based on transcript ID!
                            for c in t.CDSs.values():
                                if c.main:
                                    c.id = f"{t.id}_CDS1"
                                else:
                                    c_count += 1
                                    c.id = f"{t.id}_CDS{c_count}"
                                c.parents = [t.id]
                                cs_count = 0
                                for cs in c.CDS_segments:
                                    cs_count += 1
                                    if cs_segment_ids:
                                        cs.id = f"{c.id}_{cs_count}" 
                                    else:
                                        cs.id = c.id
                                    cs.parents = [t.id]
                                c.generate_protein()


                        if tmains > 1:
                            print(f"Error: There shouldn't be more than one main transcript for {g.id} in {self.id} annotation.")

                        e_temp = set()
                        e_parents = {}
                        for t in g.transcripts.values():
                            # gather all unique exons and sort them, then grant IDs
                            # based on start and end coordinates and start renaming
                            for e in t.exons:
                                e_unique = f"{e.start}_{e.end}_{e.strand}"
                                e_temp.add((e.start, e.end, e.strand))
                                if e_unique in e_parents:
                                    e_parents[e_unique].add(t.id)
                                else:
                                    e_parents[e_unique] = set([t.id])
                        e_temp = list(e_temp)
                        e_temp = sorted(e_temp, key=lambda x: (x[0], x[1]))
                        e_ids = {}
                        e_count = 0
                        for e in e_temp:
                            e_count += 1
                            e_count_s = "{:03d}".format(e_count)
                            e_ids[f"{e[0]}_{e[1]}_{e[2]}"] = f"{g_prefix}_e{e_count_s}"
                                
                        for t in g.transcripts.values():
                            for e in t.exons:
                                new_parents = e_parents[f"{e.start}_{e.end}_{e.strand}"]
                                e.parents = list(new_parents.copy())
                                e.parents.sort()
                                e.id = e_ids[f"{e.start}_{e.end}_{e.strand}"]

                        u_temp = set()
                        u_parents = {}
                        for t in g.transcripts.values():
                            # gather all unique UTR segments and sort them, then grant IDs
                            # based on start and end coordinates and start renaming
                            for c in t.CDSs.values():
                                for u in c.UTRs:
                                    u_unique = f"{u.start}_{u.end}"
                                    u_temp.add((u.start, u.end))
                                    if u_unique in u_parents:
                                        u_parents[u_unique].add(t.id)
                                    else:
                                        u_parents[u_unique] = set([t.id])
                        u_temp = list(u_temp)
                        u_temp = sorted(u_temp, key=lambda x: (x[0], x[1]))
                        u_ids = {}
                        u_count = 0
                        for u in u_temp:
                            u_count += 1
                            u_count_s = "{:03d}".format(u_count)
                            u_ids[f"{u[0]}_{u[1]}"] = f"{g_prefix}_u{u_count_s}"
                                
                        for t in g.transcripts.values():
                            for c in t.CDSs.values():
                                for u in c.UTRs:
                                    new_parents = u_parents[f"{u.start}_{u.end}"]
                                    u.parents = list(new_parents.copy())
                                    u.parents.sort()
                                    changed.add("UTRs")
                                    u.id = u_ids[f"{u.start}_{u.end}"]


            else:
                if correspondences:
                    correspondences_d = {}
                changed = set(["genes", "transcripts", "exons", "CDSs", "CDS_segments"])
                for genes in self.chrs.values():
                    g_count = 0
                    for g in genes.values():
                        progress_bar.update(1)
                        g_count += spacer
                        g_count_s = "{:05d}".format(g_count)
                        if suffix:
                            g_prefix = f"{prefix}{g.ch}g{g_count_s}_{suffix}"
                        else:
                            g_prefix = f"{prefix}{g.ch}g{g_count_s}"
                        if correspondences:
                            correspondences_d[g_prefix] = g.id
                        g.id = g_prefix
                        t_count = 1
                        tmains = 0
                        for t in g.transcripts.values():
                            t.parents = [g_prefix]
                            if t.main:
                                tmains += 1
                                t.id = f"{g_prefix}_t001"
                            else:
                                t_count += 1
                                t_count_s = "{:03d}".format(t_count)
                                t.id = f"{g_prefix}_t{t_count_s}"
                            c_count = 1
                            # CDSs renamed based on transcript ID!
                            for c in t.CDSs.values():
                                if c.main:
                                    c.id = f"{t.id}_CDS1"
                                else:
                                    c_count += 1
                                    c.id = f"{t.id}_CDS{c_count}"
                                c.parents = [t.id]
                                cs_count = 0
                                for cs in c.CDS_segments:
                                    cs_count += 1
                                    if cs_segment_ids:
                                        cs.id = f"{c.id}_{cs_count}" 
                                    else:
                                        cs.id = c.id
                                    cs.parents = [t.id]
                                c.generate_protein()
                        
                        if tmains > 1:
                            print(f"Error: There shouldn't be more than one main transcript for {g.id} in {self.id} annotation.")

                        e_temp = set()
                        e_parents = {}
                        for t in g.transcripts.values():
                            # gather all unique exons and sort them, then grant IDs
                            # based on start and end coordinates and start renaming
                            for e in t.exons:
                                e_unique = f"{e.start}_{e.end}_{e.strand}"
                                e_temp.add((e.start, e.end, e.strand))
                                if e_unique in e_parents:
                                    e_parents[e_unique].add(t.id)
                                else:
                                    e_parents[e_unique] = set([t.id])
                        e_temp = list(e_temp)
                        e_temp = sorted(e_temp, key=lambda x: (x[0], x[1]))
                        e_ids = {}
                        e_count = 0
                        for e in e_temp:
                            e_count += 1
                            e_count_s = "{:03d}".format(e_count)
                            e_ids[f"{e[0]}_{e[1]}_{e[2]}"] = f"{g_prefix}_e{e_count_s}"
                                
                        for t in g.transcripts.values():
                            for e in t.exons:
                                new_parents = e_parents[f"{e.start}_{e.end}_{e.strand}"]
                                e.parents = list(new_parents.copy())
                                e.parents.sort()
                                e.id = e_ids[f"{e.start}_{e.end}_{e.strand}"]

                        u_temp = set()
                        u_parents = {}
                        for t in g.transcripts.values():
                            # gather all unique UTR segments and sort them, then grant IDs
                            # based on start and end coordinates and start renaming
                            for c in t.CDSs.values():
                                for u in c.UTRs:
                                    u_unique = f"{u.start}_{u.end}"
                                    u_temp.add((u.start, u.end))
                                    if u_unique in u_parents:
                                        u_parents[u_unique].add(t.id)
                                    else:
                                        u_parents[u_unique] = set([t.id])
                        u_temp = list(u_temp)
                        u_temp = sorted(u_temp, key=lambda x: (x[0], x[1]))
                        u_ids = {}
                        u_count = 0
                        for u in u_temp:
                            u_count += 1
                            u_count_s = "{:03d}".format(u_count)
                            u_ids[f"{u[0]}_{u[1]}"] = f"{g_prefix}_u{u_count_s}"
                                
                        for t in g.transcripts.values():
                            for c in t.CDSs.values():
                                for u in c.UTRs:
                                    new_parents = u_parents[f"{u.start}_{u.end}"]
                                    u.parents = list(new_parents.copy())
                                    u.parents.sort()
                                    changed.add("UTRs")
                                    u.id = u_ids[f"{u.start}_{u.end}"]

        else:
            # Fix to change UTR names
            if not prefix:
                for genes in self.chrs.values():
                    for g in genes.values():
                        progress_bar.update(1)

                        if g.id.endswith("_gene"):
                            g_prefix = g.id[:-5]
                        elif g.id.endswith("gene"):
                            g_prefix = g.id[:-4]
                        elif g.id.startswith("gene:"):
                            g_prefix = g.id[5:]
                        elif g.id.startswith("gene-"):
                            g_prefix = g.id[5:]    
                        elif g.id.startswith("gene"):
                            g_prefix = g.id[4:]
                        else:
                            g_prefix = g.id

                        g_id_contained = False
                        g_id_not_contained = False

                        t_count = 1
                        tmains = 0
                        for t in g.transcripts.values():
                            t.parents = [g.id]
                            if g_prefix in t.id and g_prefix != t.id:
                                g_id_contained = True
                                continue
                            g_id_not_contained = True
                            if t.main:
                                tmains += 1
                                changed.add("transcripts")
                                t.id = f"{g_prefix}_t001"
                            else:
                                t_count += 1
                                t_count_s = "{:03d}".format(t_count)
                                t.id = f"{g_prefix}_t{t_count_s}"

                        if tmains > 1:
                            print(f"Error: There shouldn't be more than one main transcript for {g.id} in {self.id} annotation.")


                        if g_id_contained and g_id_not_contained:
                            print(f"Warning: {self.id} annotation has a mix of transcript formats and renaming errors could occur!")
                            
                        g_id_contained = False
                        g_id_not_contained = False
                        for t in g.transcripts.values():
                            # CDSs renamed based on transcript ID!
                            c_count = 1
                            for c in t.CDSs.values():
                                c.parents = [t.id]
                                for cs in c.CDS_segments:
                                    cs.parents = [t.id]

                                if g_prefix in c.id:
                                    g_id_contained = True
                                    continue
                                g_id_not_contained = True
                                if c.main:
                                    changed.add("CDSs")
                                    c.id = f"{t.id}_CDS1"
                                else:
                                    c_count += 1
                                    c.id = f"{t.id}_CDS{c_count}"
                                cs_count = 0
                                for cs in c.CDS_segments:
                                    changed.add("CDS_segments")
                                    cs_count += 1
                                    if cs_segment_ids:
                                        cs.id = f"{c.id}_{cs_count}" 
                                    else:
                                        cs.id = c.id
                                c.generate_protein()

                        if g_id_contained and g_id_not_contained:
                            print(f"Warning: {self.id} annotation has a mix of CDS formats and renaming errors could occur!")

                        
                        e_temp = set()
                        e_parents = {}
                        for t in g.transcripts.values():
                            # gather all unique exons and sort them, then grant IDs
                            # based on start and end coordinates and start renaming
                            for e in t.exons:
                                e_unique = f"{e.start}_{e.end}_{e.strand}"
                                e_temp.add((e.start, e.end, e.strand))
                                if e_unique in e_parents:
                                    e_parents[e_unique].add(t.id)
                                else:
                                    e_parents[e_unique] = set([t.id])
                        e_temp = list(e_temp)
                        e_temp = sorted(e_temp, key=lambda x: (x[0], x[1]))
                        e_ids = {}
                        e_count = 0
                        for e in e_temp:
                            e_count += 1
                            e_count_s = "{:03d}".format(e_count)
                            e_ids[f"{e[0]}_{e[1]}_{e[2]}"] = f"{g_prefix}_e{e_count_s}"

                        for t in g.transcripts.values():
                            for e in t.exons:
                                new_parents = e_parents[f"{e.start}_{e.end}_{e.strand}"]
                                e.parents = list(new_parents.copy())
                                e.parents.sort()
                                if g_prefix in e.id:
                                    continue
                                changed.add("exons")
                                e.id = e_ids[f"{e.start}_{e.end}_{e.strand}"]

                        u_temp = set()
                        u_parents = {}
                        for t in g.transcripts.values():
                            # gather all unique UTR segments and sort them, then grant IDs
                            # based on start and end coordinates and start renaming
                            for c in t.CDSs.values():
                                for u in c.UTRs:
                                    u_unique = f"{u.start}_{u.end}"
                                    u_temp.add((u.start, u.end))
                                    if u_unique in u_parents:
                                        u_parents[u_unique].add(t.id)
                                    else:
                                        u_parents[u_unique] = set([t.id])
                        u_temp = list(u_temp)
                        u_temp = sorted(u_temp, key=lambda x: (x[0], x[1]))
                        u_ids = {}
                        u_count = 0
                        for u in u_temp:
                            u_count += 1
                            u_count_s = "{:03d}".format(u_count)
                            u_ids[f"{u[0]}_{u[1]}"] = f"{g_prefix}_u{u_count_s}"
                                
                        for t in g.transcripts.values():
                            for c in t.CDSs.values():
                                for u in c.UTRs:
                                    new_parents = u_parents[f"{u.start}_{u.end}"]
                                    u.parents = list(new_parents.copy())
                                    u.parents.sort()
                                    if g_prefix in u.id:
                                        continue
                                    changed.add("UTRs")
                                    u.id = u_ids[f"{u.start}_{u.end}"]

            else:
                print("Error: you cannot provide a prefix and ask for minimal gene id renaming")

        progress_bar.close()
        self.update_attributes(extra_attributes=extra_attributes)
        self.update_keys()
        self.update_gene_and_transcript_list()


        if not minimal:
            self.renamed_ids = True

        if correspondences:
            out_text = "old_gene_id\tnew_gene_id\n"
            for k, v in correspondences_d.items():
                out_text += f"{v}\t{k}\n"
            f_out = open(f"{export_folder}{self.id}{self.feature_suffix}_rename_eqs.tsv", "w", encoding="utf-8")
            f_out.write(out_text)
            f_out.close()

        now = time.time()
        lapse = now - start
        print(f"\nRenaming {self.id} ids in minimal='{minimal}' and prefix='{prefix}' mode changing={list(changed)} features took {round(lapse/60, 1)} minutes")

    def update_keys(self):
        # Check if stdout or stderr are redirected to files
        stdout_redirected = not sys.stdout.isatty()
        stderr_redirected = not sys.stderr.isatty()

        # Disable tqdm if stdout or stderr are redirected
        if stdout_redirected or stderr_redirected:
            disable = True
        else:
            disable = False
        progress_bar = tqdm(total=len(self.all_gene_ids.keys())*3, disable=disable,
                                bar_format=(
                    f'\033[1;95mUpdating {self.id} dictionary keys:\033[0m '
                    '{percentage:3.0f}%|'
                    f'\033[1;95m{{bar}}\033[0m| '
                    '{n}/{total} [{elapsed}<{remaining}]'))


        # fixing CDS keys
        for genes in self.chrs.values():
            for g in genes.values():
                progress_bar.update(1)
                for t in g.transcripts.values():
                    new_CDSs = {}
                    for c in t.CDSs.values():
                        if c.id not in new_CDSs:
                            new_CDSs[c.id] = c.copy()
                        else:
                            print(f"Error: Repeated {c.id} CDS in {t.id} transcript in {self.id} annotation when updating keys.")
                    t.CDSs = new_CDSs.copy()

        for genes in self.chrs.values():
            for g in genes.values():
                progress_bar.update(1)
                new_transcripts = {}
                for t in g.transcripts.values():
                    if t.id not in new_transcripts:
                        new_transcripts[t.id] = t.copy()
                    else:
                        print(f"Error: Repeated {t.id} transcript in {g.id} gene in {self.id} annotation when updating keys.")
                g.transcripts = new_transcripts.copy()


        for chrom, genes in self.chrs.items():
            new_genes = {}
            for g in genes.values():
                progress_bar.update(1)
                new_genes[g.id] = g.copy()
            self.chrs[chrom] = new_genes.copy()
            del new_genes

        progress_bar.close()

    def update_attributes(self, clean:bool=False, featurecountsID:bool=False, aliases:bool=True, extra_attributes:bool=False, symbols:bool=False, symbols_as_descriptors=False):

        # Check if stdout or stderr are redirected to files
        stdout_redirected = not sys.stdout.isatty()
        stderr_redirected = not sys.stderr.isatty()

        # Disable tqdm if stdout or stderr are redirected
        if stdout_redirected or stderr_redirected:
            disable = True
        else:
            disable = False
        total = len(self.all_gene_ids.keys())

        if extra_attributes:
            total *= 2

        if featurecountsID:
            self.featurecounts = True

        progress_bar = tqdm(total=total, disable=disable,
                                bar_format=(
                    f'\033[1;95mUpdating {self.id} feature attributes:\033[0m '
                    '{percentage:3.0f}%|'
                    f'\033[1;95m{{bar}}\033[0m| '
                    '{n}/{total} [{elapsed}<{remaining}]'))
        for genes in self.chrs.values():
            for g in genes.values():
                progress_bar.update(1)
                g.attributes = f"ID={g.id}"
                g.attributes += f";Name={g.id}"
                new_symbols = ",".join(g.symbols)
                if new_symbols != "" and symbols:
                    g.attributes += f";Symbol={new_symbols}"
                if new_symbols != "" and symbols_as_descriptors:
                    g.attributes += f";Description={new_symbols}"
                if not clean and g.misc_attributes != "":
                    g.attributes += f";{g.misc_attributes}"
                g.attributes += f";pseudogene={g.pseudogene};transposable={g.transposable}"
                if featurecountsID:
                    g.attributes += f";featurecounts_id={g.id}"
                if aliases and g.aliases != []:
                    new_aliases = ",".join(g.aliases)
                    g.attributes += f";Alias={new_aliases}"
                for t in g.transcripts.values():
                    new_parents = ",".join(t.parents)
                    t.attributes = f"ID={t.id};Parent={new_parents}"
                    if not clean and t.misc_attributes != "":
                        t.attributes += f";{t.misc_attributes}"
                    if featurecountsID:
                        t.attributes += f";featurecounts_id={g.id}"
                    if aliases and t.aliases != []:
                        new_aliases = ",".join(t.aliases)
                        t.attributes += f";Alias={new_aliases}"
                    for e in t.exons:
                        new_parents = ",".join(e.parents)
                        e.attributes = f"ID={e.id};Parent={new_parents}"
                        if not clean and e.misc_attributes != "":
                            e.attributes += f";{e.misc_attributes}"
                        if featurecountsID:
                            e.attributes += f";featurecounts_id={g.id}"
                    for c in t.CDSs.values():
                        new_parents = ",".join(c.parents)
                        c.attributes = f"ID={c.id};Parent={new_parents}"
                        if not clean and c.misc_attributes != "":
                            c.attributes += f";{c.misc_attributes}"
                        if featurecountsID:
                            c.attributes += f";featurecounts_id={g.id}"
                        for cs in c.CDS_segments:
                            new_parents = ",".join(cs.parents)
                            cs.attributes = f"ID={cs.id};Parent={new_parents}"
                            if not clean and cs.misc_attributes != "":
                                cs.attributes += f";{cs.misc_attributes}"
                            if featurecountsID:
                                cs.attributes += f";featurecounts_id={g.id}"
                        for u in c.UTRs:
                            new_parents = ",".join(u.parents)
                            u.attributes = f"ID={u.id};Parent={new_parents}"
                            if not clean and u.misc_attributes != "":
                                u.attributes += f";{u.misc_attributes}"
                            if featurecountsID:
                                u.attributes += f";featurecounts_id={g.id}"
        
        if extra_attributes:
            for genes in self.chrs.values():
                for g in genes.values():
                    progress_bar.update(1)
                    
                    g.attributes += f";reliable_score={g.reliable_score}"
                    g.attributes += f";remove={g.remove}"
                    g.attributes += f";rescue={g.rescue}"
                    blasts = []
                    for b in g.blast_hits:
                        blasts.append(f"{b.source}_{b.score}")
                    blasts = ",".join(blasts)
                    if blasts != "":
                        g.attributes += f";blasts={blasts}"
                    alternative_transcript_rescue = ",".join(list(g.alternative_transcript_rescue))
                    if alternative_transcript_rescue != "":
                        g.attributes += f";alternative_transcript_rescue={alternative_transcript_rescue}"
                    overlaps = []
                    for o in g.overlaps["self"]:
                        if o.score >= 5:
                            overlaps.append(o.id)
                    overlaps = ",".join(overlaps)
                    if overlaps != "":
                        g.attributes += f";CDS_orientated_overlaps={overlaps}"

                    gene_masked_fraction = g.masked_fraction
                    transcript_masked_fraction = 0
                    CDS_masked_fraction = 0
                    gene_GC_content = g.gc_content
                    transcript_GC_content = 0
                    CDS_GC_content = 0

                    for t in g.transcripts.values():
                        if t.main:
                            transcript_masked_fraction = t.masked_fraction
                            transcript_GC_content = t.gc_content
                            for c in t.CDSs.values():
                                if c.main:
                                    CDS_masked_fraction = c.masked_fraction
                                    CDS_GC_content = c.gc_content

                    g.attributes += f";gene_masked_fraction={gene_masked_fraction}"
                    g.attributes += f";transcript_masked_fraction={transcript_masked_fraction}"
                    g.attributes += f";CDS_masked_fraction={CDS_masked_fraction}"
                    g.attributes += f";gene_GC_content={gene_GC_content}"
                    g.attributes += f";transcript_GC_content={transcript_GC_content}"
                    g.attributes += f";CDS_GC_content={CDS_GC_content}"
                    g.attributes += f";intron_nested={g.intron_nested}"
                    g.attributes += f";intron_nested_fully_contained={g.intron_nested_fully_contained}"
                    g.attributes += f";intron_nested_single={g.intron_nested_single}"
                    g.attributes += f";intron_UTR_nested={g.UTR_intron_nested}"

        progress_bar.close()

        if self.renamed_ids:
            self.aegis = True

        if clean:
            self.clean = True

        self.update_suffixes()

    def temporary_add_missing_properties(self):
        pass

    def add_blast_hits(self, source, blastfile, mode:str="protein"):
        """
        Blast results could be at a gene, transcript or protein level. In future,
        other features can also be easily included.

        This function for now assumes one hit per protein!
        """
        accepted_levels = ["protein", "CDS", "transcript", "gene"]

        if mode in accepted_levels:
            hits = {}
            f_in = open(blastfile, encoding="utf-8")
            for line in f_in:
                line = line.strip().split("\t")
                ID = line[0]
                evalue = float(line[-2])
                bitscore = float(line[-1])
                if ID not in hits:
                    hits[ID] = BlastHit(source, bitscore, evalue)
                else:
                    print(f"Warning: for now the add blast hits method does not accept more than one blast hit per ID ({ID})")
            f_in.close()
            if mode == "protein":
                for protein_id in self.protein_equivalences:
                    chr, g, t, c = self.all_protein_ids[protein_id]
                    if protein_id in hits:
                        self.chrs[chr][g].transcripts[t].CDSs[c].protein.blast_hits.append(hits[protein_id])
                    for protein_id_copy in self.protein_equivalences[protein_id]:
                        chr, g, t, c = self.all_protein_ids[protein_id_copy]
                        if protein_id_copy in self.chrs[chr] and protein_id in hits:
                            self.chrs[chr][g].transcripts[t].CDSs[c].protein.blast_hits.append(hits[protein_id])
            else:
                "Adding blast hits is not available for gene, transcripts and CDSs yet."
        else:
            print(f"Warning: {mode} chosen is not in accepted list of choices=['protein', 'CDS', 'transcript', 'gene]")

    def add_qualitative_info_to_overlaps(self):
        """
        Number of unique full segment overlaps between genes including all transcript variants.
        """
        # Check if stdout or stderr are redirected to files
        stdout_redirected = not sys.stdout.isatty()
        stderr_redirected = not sys.stderr.isatty()

        # Disable tqdm if stdout or stderr are redirected
        if stdout_redirected or stderr_redirected:
            disable = True
        else:
            disable = False
        total = len(self.all_gene_ids.keys())

        progress_bar = tqdm(total=total, disable=disable,
                                bar_format=(
                    f'\033[1;95mAdding qualitative info to {self.id} overlaps:\033[0m '
                    '{percentage:3.0f}%|'
                    f'\033[1;95m{{bar}}\033[0m| '
                    '{n}/{total} [{elapsed}<{remaining}]'))
        for genes in self.chrs.values():
            for g in genes.values():
                progress_bar.update(1)
                for o in g.overlaps["self"]:
                    if o.exon_query_percent > 0:

                        g_exons, g_CDSs, g_UTRs = set(), set(), set()
                        o_exons, o_CDSs, o_UTRs = set(), set(), set()

                        for t in g.transcripts.values():
                            if t.main:
                                g_exons.update(set([(e.start, e.end) for e in t.exons]))
                                for c in t.CDSs.values():
                                    if c.main:
                                        g_CDSs.update(set([(cs.start, cs.end, cs.frame) for cs in c.CDS_segments]))
                                        g_UTRs.update(set([(u.start, u.end) for u in c.UTRs]))
                                
                        chrom = self.all_gene_ids[o.id]
                        for t in self.chrs[chrom][o.id].transcripts.values():
                            if t.main:
                                o_exons.update(set([(e.start, e.end) for e in t.exons]))
                                for c in t.CDSs.values():
                                    if c.main:
                                        o_CDSs.update(set([(cs.start, cs.end, cs.frame) for cs in c.CDS_segments]))
                                        o_UTRs.update(set([(u.start, u.end) for u in c.UTRs]))

                        for e1 in g_exons:
                            for e2 in o_exons:
                                if e1 == e2:
                                    o.full_exon_overlaps += 1

                        for c1 in g_CDSs:
                            for c2 in o_CDSs:
                                if c1[0] == c2[0] and c1[1] == c2[1]:
                                    o.full_CDS_overlaps += 1
                                    if c1[2] == c2[2]:
                                        o.full_protein_overlaps += 1

                        for u1 in g_UTRs:
                            for u2 in o_UTRs:
                                if u1 == u2:
                                    o.full_UTR_overlaps += 1

        progress_bar.close()

    def add_qualitative_info_to_overlaps_old(self):

        # Check if stdout or stderr are redirected to files
        stdout_redirected = not sys.stdout.isatty()
        stderr_redirected = not sys.stderr.isatty()

        # Disable tqdm if stdout or stderr are redirected
        if stdout_redirected or stderr_redirected:
            disable = True
        else:
            disable = False
        total = len(self.all_gene_ids.keys())

        progress_bar = tqdm(total=total, disable=disable,
                                bar_format=(
                    f'\033[1;95mAdding qualitative info to {self.id} overlaps:\033[0m '
                    '{percentage:3.0f}%|'
                    f'\033[1;95m{{bar}}\033[0m| '
                    '{n}/{total} [{elapsed}<{remaining}]'))
        for genes in self.chrs.values():
            for g in genes.values():
                progress_bar.update(1)
                for o in g.overlaps["self"]:
                    if o.exon_query_percent > 0:
                        g_exons = []
                        g_CDSs = []
                        g_UTRs = []
                        for t in g.transcripts.values():
                            if t.main:
                                g_exons = t.exons.copy()
                                for c in t.CDSs.values():
                                    if c.main:
                                        g_CDSs = c.CDS_segments.copy()
                                        g_UTRs = c.UTRs.copy()
                        chrom = self.all_gene_ids[o.id]
                        o_exons = []
                        o_CDSs = []
                        o_UTRs = []
                        for t in self.chrs[chrom][o.id].transcripts.values():
                            if t.main:
                                o_exons = t.exons.copy()
                                for c in t.CDSs.values():
                                    if c.main:
                                        o_CDSs = c.CDS_segments.copy()
                                        o_UTRs = c.UTRs.copy()

                        for e1 in g_exons:
                            for e2 in o_exons:
                                if e1.equal_coordinates(e2):
                                    o.full_exon_overlaps += 1

                        for c1 in g_CDSs:
                            for c2 in o_CDSs:
                                if c1.equal_coordinates(c2):
                                    o.full_CDS_overlaps += 1

                        for u1 in g_UTRs:
                            for u2 in o_UTRs:
                                if u1.equal_coordinates(u2):
                                    o.full_UTR_overlaps += 1
        progress_bar.close()

    def clear_overlaps_with_selected_CDSs(self):
        for genes in self.chrs.values():
            for g in genes.values():
                g.overlap_with_selected_CDS = False   

    def clear_overlaps_with_selected_exons(self):
        for genes in self.chrs.values():
            for g in genes.values():
                g.overlap_with_selected_exon = False        

    def mark_intron_nesting(self, ignore_removed:bool=True):
        for genes in self.chrs.values():
            for g in genes.values():
                for o in g.overlaps["self"]:
                    if o.score < 5:
                        chrom = self.all_gene_ids[o.id]
                        if ignore_removed:
                            if self.chrs[chrom][o.id].remove and not self.chrs[chrom][o.id].rescue:
                                continue
                        if (o.exon_query_percent == 0 and o.exon_target_percent == 0) and (g.start > self.chrs[chrom][o.id].start or g.end < self.chrs[chrom][o.id].end):
                            g.intron_nested = True
                            if self.chrs[chrom][o.id].start < g.start and self.chrs[chrom][o.id].end > g.end:
                                g.intron_nested_fully_contained = True
                            
                            c_start_target = None
                            c_end_target = None
                            c_start_query = None
                            c_end_query = None

                            for t in self.chrs[chrom][o.id].transcripts.values():
                                if t.main:
                                    for c in t.CDSs.values():
                                        if c.main:
                                            c_start_target = c.start
                                            c_end_target = c.end

                            for t in g.transcripts.values():
                                if t.main:
                                    for c in t.CDSs.values():
                                        if c.main:
                                            c_start_query = c.start
                                            c_end_query = c.end

                            # UTR intron nested means that a main CDS of a gene finishes and starts outside of the overlaped gene's CDS region
                            if c_start_target != None and c_start_query != None:
                                if c_start_query > c_end_target or c_end_query < c_start_target:
                                    g.UTR_intron_nested = True
                                    
                            for t in self.chrs[chrom][o.id].transcripts.values():
                                if t.main:
                                    for i in t.introns:
                                        if i.start < g.start and i.end > g.end:
                                            g.intron_nested_single = True

    def mark_noisy_genes(self, protein_size:int=50, intron_size:int=100000, remove_noncoding:bool=True, remove_masked:bool=True):
        # Check if stdout or stderr are redirected to files
        stdout_redirected = not sys.stdout.isatty()
        stderr_redirected = not sys.stderr.isatty()

        # Disable tqdm if stdout or stderr are redirected
        if stdout_redirected or stderr_redirected:
            disable = True
        else:
            disable = False
        progress_bar = tqdm(total=len(self.all_gene_ids.keys()), disable=disable,
                                bar_format=(
                    f'\033[1;91mMark {self.id} noisy genes:\033[0m '
                    '{percentage:3.0f}%|'
                    f'\033[1;91m{{bar}}\033[0m| '
                    '{n}/{total} [{elapsed}<{remaining}]'))
        for genes in self.chrs.values():
            for g in genes.values():
                progress_bar.update(1)
                if remove_masked:
                    if g.masked_fraction == 1:
                        g.remove = True
                if remove_noncoding:
                    if not g.coding:
                        g.remove = True
                    else:
                        for t in g.transcripts.values():
                            if t.main:
                                if t.masked_fraction == 1:
                                    g.remove = True
                                for c in t.CDSs.values():
                                    if c.main:
                                        if c.size < (protein_size * 3) or c.masked_fraction == 1:
                                            g.remove = True
                for t in g.transcripts.values():
                    if t.main:
                        for i in t.introns:
                            if i.size > intron_size:
                                g.remove = True
        progress_bar.close()

    def add_reliable_CDS_evidence_score(self):
        # Check if stdout or stderr are redirected to files
        stdout_redirected = not sys.stdout.isatty()
        stderr_redirected = not sys.stderr.isatty()

        # Disable tqdm if stdout or stderr are redirected
        if stdout_redirected or stderr_redirected:
            disable = True
        else:
            disable = False
        progress_bar = tqdm(total=len(self.all_gene_ids.keys()), disable=disable,
                                bar_format=(
                    f'\033[1;91mMark {self.id} reliable CDS evidences:\033[0m '
                    '{percentage:3.0f}%|'
                    f'\033[1;91m{{bar}}\033[0m| '
                    '{n}/{total} [{elapsed}<{remaining}]'))

        for genes in self.chrs.values():
            for g in genes.values():
                progress_bar.update(1)
                if not g.remove:
                    for o in g.overlaps["self"]:
                        if not self.chrs[g.ch][o.id].remove:
                            if o.score == 11:
                                g.reliable_score += 1

        for genes in self.chrs.values():
            for g in genes.values():
                if g.reliable_score == 0:
                    g.remove = True
        progress_bar.close()

    def mark_reliable_CDS_evidences(self, unreliable_sources:list=["GlimmerHMM", "geneid_v1.4"]):
        # Check if stdout or stderr are redirected to files
        stdout_redirected = not sys.stdout.isatty()
        stderr_redirected = not sys.stderr.isatty()

        # Disable tqdm if stdout or stderr are redirected
        if stdout_redirected or stderr_redirected:
            disable = True
        else:
            disable = False
        progress_bar = tqdm(total=len(self.all_gene_ids.keys()), disable=disable,
                                bar_format=(
                    f'\033[1;91mMark {self.id} reliable CDS evidences:\033[0m '
                    '{percentage:3.0f}%|'
                    f'\033[1;91m{{bar}}\033[0m| '
                    '{n}/{total} [{elapsed}<{remaining}]'))

        for genes in self.chrs.values():
            for g in genes.values():
                progress_bar.update(1)
                if not g.remove:
                    for o in g.overlaps["self"]:
                        if not self.chrs[g.ch][o.id].remove:
                            if (self.chrs[g.ch][o.id].source not in unreliable_sources) or (g.source not in unreliable_sources):
                                if o.score == 11:
                                    g.reliable = True
                                    if not self.chrs[g.ch][o.id].reliable:
                                        self.chrs[g.ch][o.id].remove = True
                                    self.chrs[g.ch][o.id].reliable = True

        for genes in self.chrs.values():
            for g in genes.values():
                if not g.reliable:
                    g.remove = True
        progress_bar.close()

    def remove_genes(self):
        # Check if stdout or stderr are redirected to files
        stdout_redirected = not sys.stdout.isatty()
        stderr_redirected = not sys.stderr.isatty()

        # Disable tqdm if stdout or stderr are redirected
        if stdout_redirected or stderr_redirected:
            disable = True
        else:
            disable = False
        progress_bar = tqdm(total=len(self.all_gene_ids.keys()), disable=disable,
                                bar_format=(
                    f'\033[1;91mRemoving {self.id} genes:\033[0m '
                    '{percentage:3.0f}%|'
                    f'\033[1;91m{{bar}}\033[0m| '
                    '{n}/{total} [{elapsed}<{remaining}]'))

        genes_to_remove = set()
        for genes in self.chrs.values():
            for g in genes.values():
                progress_bar.update(1)
                if g.remove and not g.rescue:
                    genes_to_remove.add(g.id)

        for g_id in genes_to_remove:
            chrom = self.all_gene_ids[g_id]
            del self.chrs[chrom][g_id]
            del self.all_gene_ids[g_id]
        progress_bar.close()
        self.remove_missing_genes_in_overlaps()

    def remove_missing_genes_in_overlaps(self):
        # Check if stdout or stderr are redirected to files
        stdout_redirected = not sys.stdout.isatty()
        stderr_redirected = not sys.stderr.isatty()

        # Disable tqdm if stdout or stderr are redirected
        if stdout_redirected or stderr_redirected:
            disable = True
        else:
            disable = False
        progress_bar = tqdm(total=len(self.all_gene_ids.keys()), disable=disable,
                                bar_format=(
                    f'\033[1;91mRemoving {self.id} missing genes in overlaps:\033[0m '
                    '{percentage:3.0f}%|'
                    f'\033[1;91m{{bar}}\033[0m| '
                    '{n}/{total} [{elapsed}<{remaining}]'))
        for genes in self.chrs.values():
            for g in genes.values():
                progress_bar.update(1)
                new_overlaps = []
                for o in g.overlaps["self"]:
                    if o.id in self.all_gene_ids:
                        new_overlaps.append(o)
                g.overlaps["self"] = new_overlaps
        progress_bar.close()              

    def mark_transcriptomic_supported_genes(self):
        # Check if stdout or stderr are redirected to files
        stdout_redirected = not sys.stdout.isatty()
        stderr_redirected = not sys.stderr.isatty()

        # Disable tqdm if stdout or stderr are redirected
        if stdout_redirected or stderr_redirected:
            disable = True
        else:
            disable = False
        progress_bar = tqdm(total=len(self.all_gene_ids.keys()), disable=disable,
                                bar_format=(
                    f'\033[1;91mMark {self.id} transcriptomic supported genes:\033[0m '
                    '{percentage:3.0f}%|'
                    f'\033[1;91m{{bar}}\033[0m| '
                    '{n}/{total} [{elapsed}<{remaining}]'))
        for genes in self.chrs.values():
            for g in genes.values():
                g.transcriptomic_evidence = False
                progress_bar.update(1)
                if "stringtie" in g.source or "psiclass" in g.source:
                    g.transcriptomic_evidence = True
                else:
                    for o in g.overlaps["self"]:
                        if o.score >= 5:
                            if o.CDS_query_percent > 30:
                                if "stringtie" in self.chrs[g.ch][o.id].source or "psiclass" in self.chrs[g.ch][o.id].source:
                                    g.transcriptomic_evidence = True
        progress_bar.close()

    def mark_abinitio_supported_genes(self, reliable_sources:list=["AUGUSTUS", "GeneMark.hmm3"]):
        # Check if stdout or stderr are redirected to files
        stdout_redirected = not sys.stdout.isatty()
        stderr_redirected = not sys.stderr.isatty()

        # Disable tqdm if stdout or stderr are redirected
        if stdout_redirected or stderr_redirected:
            disable = True
        else:
            disable = False
        progress_bar = tqdm(total=len(self.all_gene_ids.keys()), disable=disable,
                                bar_format=(
                    f'\033[1;91mMark {self.id} abinitio supported genes:\033[0m '
                    '{percentage:3.0f}%|'
                    f'\033[1;91m{{bar}}\033[0m| '
                    '{n}/{total} [{elapsed}<{remaining}]'))
        for genes in self.chrs.values():
            for g in genes.values():
                g.abinitio_evidence = False
                progress_bar.update(1)
                if g.source in reliable_sources:
                    g.abinitio_evidence = True
                else:
                    for o in g.overlaps["self"]:
                        if o.score >= 5:
                            if o.CDS_query_percent > 30:
                                if self.chrs[g.ch][o.id].source in reliable_sources:
                                    g.abinitio_evidence = True
        progress_bar.close()

    def mark_overlap_with_reliable_genes(self):
                # Check if stdout or stderr are redirected to files
        stdout_redirected = not sys.stdout.isatty()
        stderr_redirected = not sys.stderr.isatty()

        # Disable tqdm if stdout or stderr are redirected
        if stdout_redirected or stderr_redirected:
            disable = True
        else:
            disable = False
        progress_bar = tqdm(total=len(self.all_gene_ids.keys()), disable=disable,
                                bar_format=(
                    f'\033[1;91mMark {self.id} overlap with reliable genes:\033[0m '
                    '{percentage:3.0f}%|'
                    f'\033[1;91m{{bar}}\033[0m| '
                    '{n}/{total} [{elapsed}<{remaining}]'))
        for genes in self.chrs.values():
            for g in genes.values():
                progress_bar.update(1)
                for o in g.overlaps["self"]:
                    if o.score >= 5 or (o.score == 1 and o.antiscore >= 5):
                        if not self.chrs[g.ch][o.id].remove:
                            g.overlap_reliable = True
        progress_bar.close()

    def find_best_gene_model(self, source_priority:list, just_with_reliables:bool=True):
        # Check if stdout or stderr are redirected to files
        stdout_redirected = not sys.stdout.isatty()
        stderr_redirected = not sys.stderr.isatty()

        # Disable tqdm if stdout or stderr are redirected
        if stdout_redirected or stderr_redirected:
            disable = True
        else:
            disable = False
        progress_bar = tqdm(total=len(self.all_gene_ids.keys()), disable=disable,
                                bar_format=(
                    f'\033[1;91mMark {self.id} noisy genes:\033[0m '
                    '{percentage:3.0f}%|'
                    f'\033[1;91m{{bar}}\033[0m| '
                    '{n}/{total} [{elapsed}<{remaining}]'))
        
        if just_with_reliables:
            for genes in self.chrs.values():
                for g in genes.values():
                    progress_bar.update(1)
                    if not g.remove:
                        for o in g.overlaps["self"]:
                            if o.score >= 5 or (o.score == 1 and o.antiscore >= 5):
                                if not self.chrs[g.ch][o.id].remove and not g.remove:
                                    if g.reliable_score > self.chrs[g.ch][o.id].reliable_score:
                                        self.chrs[g.ch][o.id].remove = True
                                    elif self.chrs[g.ch][o.id].reliable_score > g.reliable_score:
                                        g.remove = True
                                    else:
                                        query_best = g.compare_protein_blast_hits(self.chrs[g.ch][o.id], source_priority)
                                        if query_best:
                                            self.chrs[g.ch][o.id].remove = True
                                        else:
                                            g.remove = True

        else:
            for genes in self.chrs.values():
                for g in genes.values():
                    progress_bar.update(1)
                    if g.remove and not g.overlap_reliable and g.transcriptomic_evidence and not g.unrescuable:
                        for o in g.overlaps["self"]:
                            if not self.chrs[g.ch][o.id].remove:
                                if o.score == 1:
                                    for t3 in g.transcripts.values():
                                        for t4 in self.chrs[g.ch][o.id].transcripts.values():
                                            for e1 in t3.exons:
                                                for c1 in t4.CDSs.values():
                                                    for cs1 in c1.CDS_segments:
                                                        overlapping, _ = overlap(e1, cs1)
                                                        if overlapping:
                                                            g.unrescuable = True
                                                            g.rescue = False
                                                            break

                            if self.chrs[g.ch][o.id].remove and not self.chrs[g.ch][o.id].overlap_reliable and self.chrs[g.ch][o.id].transcriptomic_evidence and not self.chrs[g.ch][o.id].unrescuable and not g.unrescuable:
                                if o.score >= 5 or (o.score == 1 and o.antiscore >= 5):
                                    query_best = g.compare_protein_blast_hits(self.chrs[g.ch][o.id], source_priority)
                                    if query_best:
                                        g.rescue = True
                                        self.chrs[g.ch][o.id].rescue = False
                                        self.chrs[g.ch][o.id].remove = True
                                        self.chrs[g.ch][o.id].unrescuable = True
                                    else:
                                        self.chrs[g.ch][o.id].rescue = True
                                        g.rescue = False
                                        g.remove = True
                                        g.unrescuable = True    

        progress_bar.close()

    def mark_overlap_with_other_selected_exons(self):
        self.clear_overlaps_with_selected_exons()
        # Check if stdout or stderr are redirected to files
        stdout_redirected = not sys.stdout.isatty()
        stderr_redirected = not sys.stderr.isatty()

        # Disable tqdm if stdout or stderr are redirected
        if stdout_redirected or stderr_redirected:
            disable = True
        else:
            disable = False
        progress_bar = tqdm(total=len(self.all_gene_ids.keys()), disable=disable,
                                bar_format=(
                    f'\033[1;91mMark {self.id} overlap with other selected exons:\033[0m '
                    '{percentage:3.0f}%|'
                    f'\033[1;91m{{bar}}\033[0m| '
                    '{n}/{total} [{elapsed}<{remaining}]'))   
        for genes in self.chrs.values():
            for g in genes.values():
                progress_bar.update(1)
                if not g.remove or g.rescue:
                    for o in g.overlaps["self"]:
                        if o.score < 5 and (self.chrs[g.ch][o.id].rescue or not self.chrs[g.ch][o.id].remove):
                            for t in g.transcripts.values():
                                for t2 in self.chrs[g.ch][o.id].transcripts.values():
                                    for e1 in t.exons:
                                        for e2 in t2.exons:
                                            overlapping, _ = overlap(e1, e2)
                                            if overlapping:
                                                g.overlap_with_selected_exon = True
                                                break
                    for o in g.overlaps["self"]:
                        if o.score == 11:
                            for o2 in self.chrs[g.ch][o.id].overlaps["self"]:
                                if o2.score < 5 and (not self.chrs[g.ch][o2.id].remove or self.chrs[g.ch][o2.id].rescue) and o2.id != g.id:
                                    for t in self.chrs[g.ch][o.id].transcripts.values():
                                        for t2 in self.chrs[g.ch][o2.id].transcripts.values():
                                            for e1 in t.exons:
                                                for e2 in t2.exons:
                                                    overlapping, _ = overlap(e1, e2)
                                                    if overlapping:
                                                        self.chrs[g.ch][o.id].overlap_with_selected_exon = True
                                                        break
        progress_bar.close()        

    def mark_overlap_with_other_selected_CDSs(self):
        self.clear_overlaps_with_selected_CDSs()
        # Check if stdout or stderr are redirected to files
        stdout_redirected = not sys.stdout.isatty()
        stderr_redirected = not sys.stderr.isatty()

        # Disable tqdm if stdout or stderr are redirected
        if stdout_redirected or stderr_redirected:
            disable = True
        else:
            disable = False
        progress_bar = tqdm(total=len(self.all_gene_ids.keys()), disable=disable,
                                bar_format=(
                    f'\033[1;91mMark {self.id} overlap with other selected CDSs:\033[0m '
                    '{percentage:3.0f}%|'
                    f'\033[1;91m{{bar}}\033[0m| '
                    '{n}/{total} [{elapsed}<{remaining}]'))        
        for genes in self.chrs.values():
            for g in genes.values():
                progress_bar.update(1)
                if not g.remove or g.rescue:
                    for o in g.overlaps["self"]:
                        if o.score < 5 and (self.chrs[g.ch][o.id].rescue or not self.chrs[g.ch][o.id].remove):
                            for t in g.transcripts.values():
                                for t2 in self.chrs[g.ch][o.id].transcripts.values():
                                    for e in t.exons:
                                        for c in t2.CDSs.values():
                                            for cs in c.CDS_segments:
                                                overlapping, _ = overlap(e, cs)
                                                if overlapping:
                                                    g.overlap_with_selected_CDS = True
                                                    break
                    for o in g.overlaps["self"]:
                        if o.score == 11:
                            for o2 in self.chrs[g.ch][o.id].overlaps["self"]:
                                if o2.score < 5 and (not self.chrs[g.ch][o2.id].remove or self.chrs[g.ch][o2.id].rescue) and o2.id != g.id:
                                    for t in self.chrs[g.ch][o.id].transcripts.values():
                                        for t2 in self.chrs[g.ch][o2.id].transcripts.values():
                                            for e in t.exons:
                                                for c in t2.CDSs.values():
                                                    for cs in c.CDS_segments:
                                                        overlapping, _ = overlap(e, cs)
                                                        if overlapping:
                                                            self.chrs[g.ch][o.id].overlap_with_selected_CDS = True
                                                            break
        progress_bar.close()

    def select_best_possible_non_overlapping_UTR(self, exon=False):
        # Check if stdout or stderr are redirected to files
        stdout_redirected = not sys.stdout.isatty()
        stderr_redirected = not sys.stderr.isatty()

        # Disable tqdm if stdout or stderr are redirected
        if stdout_redirected or stderr_redirected:
            disable = True
        else:
            disable = False
        progress_bar = tqdm(total=len(self.all_gene_ids.keys()), disable=disable,
                                bar_format=(
                    f'\033[1;91mSelect {self.id} best possible non overlapping UTR:\033[0m '
                    '{percentage:3.0f}%|'
                    f'\033[1;91m{{bar}}\033[0m| '
                    '{n}/{total} [{elapsed}<{remaining}]'))
        
        if not exon:
            for genes in self.chrs.values():
                for g in genes.values():
                    g.update()
                    progress_bar.update(1)
                    if not g.overlap_with_selected_CDS:
                        sizes = [g.size]
                    else:
                        sizes = []
                    if not g.remove or g.rescue:
                        if len(g.transcripts) < 2:
                            for o in g.overlaps["self"]:
                                if o.score == 11:
                                    if not self.chrs[g.ch][o.id].overlap_with_selected_CDS:
                                        sizes.append(self.chrs[g.ch][o.id].size)

                            rescue_transcripts = {}
                            if sizes != []:
                                max_gene_size = max(sizes)
                                for o in g.overlaps["self"]:
                                    if o.score == 11:
                                        if not self.chrs[g.ch][o.id].overlap_with_selected_CDS and self.chrs[g.ch][o.id].size == max_gene_size and rescue_transcripts == {}:
                                            for t in self.chrs[g.ch][o.id].transcripts.values():
                                                rescue_transcripts[t.id] = t.copy()

                            else:
                                sizes = [g.size]
                                for o in g.overlaps["self"]:
                                    if o.score == 11:
                                        sizes.append(self.chrs[g.ch][o.id].size)
                                min_gene_size = min(sizes)
                                for o in g.overlaps["self"]:
                                    if o.score == 11:
                                        if not self.chrs[g.ch][o.id].overlap_with_selected_CDS and self.chrs[g.ch][o.id].size == min_gene_size and rescue_transcripts == {}:
                                            for t in self.chrs[g.ch][o.id].transcripts.values():
                                                rescue_transcripts[t.id] = t.copy()

                            if rescue_transcripts != {}:
                                g.transcripts = rescue_transcripts.copy()
                                g.update()
        else:
            for genes in self.chrs.values():
                for g in genes.values():
                    g.update()
                    progress_bar.update(1)
                    if not g.overlap_with_selected_exon:
                        sizes = [g.size]
                    else:
                        sizes = []
                    if not g.remove or g.rescue:
                        if len(g.transcripts) < 2:
                            for o in g.overlaps["self"]:
                                if o.score == 11:
                                    if not self.chrs[g.ch][o.id].overlap_with_selected_exon:
                                        sizes.append(self.chrs[g.ch][o.id].size)

                            rescue_transcripts = {}
                            if sizes != []:
                                max_gene_size = max(sizes)
                                for o in g.overlaps["self"]:
                                    if o.score == 11:
                                        if not self.chrs[g.ch][o.id].overlap_with_selected_exon and self.chrs[g.ch][o.id].size == max_gene_size and rescue_transcripts == {}:
                                            for t in self.chrs[g.ch][o.id].transcripts.values():
                                                rescue_transcripts[t.id] = t.copy()

                            else:
                                sizes = [g.size]
                                for o in g.overlaps["self"]:
                                    if o.score == 11:
                                        sizes.append(self.chrs[g.ch][o.id].size)
                                min_gene_size = min(sizes)
                                for o in g.overlaps["self"]:
                                    if o.score == 11:
                                        if not self.chrs[g.ch][o.id].overlap_with_selected_exon and self.chrs[g.ch][o.id].size == min_gene_size and rescue_transcripts == {}:
                                            for t in self.chrs[g.ch][o.id].transcripts.values():
                                                rescue_transcripts[t.id] = t.copy()

                            if rescue_transcripts != {}:
                                g.transcripts = rescue_transcripts.copy()      
                                g.update()      

        progress_bar.close()
        self.rename_ids()
        self.remove_duplicate_transcripts()
        self.update()

    def remove_duplicate_transcripts(self):
        # Check if stdout or stderr are redirected to files
        stdout_redirected = not sys.stdout.isatty()
        stderr_redirected = not sys.stderr.isatty()

        # Disable tqdm if stdout or stderr are redirected
        if stdout_redirected or stderr_redirected:
            disable = True
        else:
            disable = False
        progress_bar = tqdm(total=len(self.all_gene_ids.keys()), disable=disable,
                                bar_format=(
                    f'\033[1;91mRemoving repeat transcripts per gene of {self.id}:\033[0m '
                    '{percentage:3.0f}%|'
                    f'\033[1;91m{{bar}}\033[0m| '
                    '{n}/{total} [{elapsed}<{remaining}]'))
        for genes in self.chrs.values():
            for g in genes.values():
                progress_bar.update(1)
                to_eliminate = []
                for t1 in g.transcripts.values():
                    if t1.id not in to_eliminate:
                        for t2 in g.transcripts.values():
                            if t2.id not in to_eliminate:
                                if t1.id == t2.id:
                                    continue
                                if t1.almost_equal(t2):
                                    to_eliminate.append(t1.id)
                to_eliminate = set(to_eliminate)
                for t_eliminate in to_eliminate:
                    del g.transcripts[t_eliminate]
        progress_bar.close()
        self.update()

    def add_better_ab_initio_models_as_alternative_transcripts(self, source_priority, reliable_sources:list=["AUGUSTUS", "GeneMark.hmm3"]):
        # Check if stdout or stderr are redirected to files
        stdout_redirected = not sys.stdout.isatty()
        stderr_redirected = not sys.stderr.isatty()

        # Disable tqdm if stdout or stderr are redirected
        if stdout_redirected or stderr_redirected:
            disable = True
        else:
            disable = False
        progress_bar = tqdm(total=len(self.all_gene_ids.keys())*3, disable=disable,
                                bar_format=(
                    f'\033[1;91mSelecting {self.id} alternative transcripts:\033[0m '
                    '{percentage:3.0f}%|'
                    f'\033[1;91m{{bar}}\033[0m| '
                    '{n}/{total} [{elapsed}<{remaining}]'))

        for genes in self.chrs.values():
            for g in genes.values():
                progress_bar.update(1)
                g_coding_ratio = 0
                full_UTR_exons = 0
                for t in g.transcripts.values():
                    if t.main:
                        g_coding_ratio = t.coding_ratio
                        for c in t.CDSs.values():
                            if c.main:
                                full_UTR_exons = c.full_UTR_exons

                if (not g.remove or g.rescue) and g_coding_ratio < 0.7 and full_UTR_exons > 0:
                    for o in g.overlaps["self"]:
                        if o.score >= 5:
                            if self.chrs[g.ch][o.id].remove and not self.chrs[g.ch][o.id].rescue and self.chrs[g.ch][o.id].source in reliable_sources:
                                query_best = g.compare_protein_blast_hits(self.chrs[g.ch][o.id], source_priority)
                                if not query_best:
                                    g.alternative_transcript_rescue.add(o.id)
        gene_groups = []

        for genes in self.chrs.values():
            for g in genes.values():
                progress_bar.update(1)
                if (not g.remove or g.rescue) and g.alternative_transcript_rescue != set():
                    for g2 in genes.values():
                        if g.id == g2.id:
                            continue
                        if (not g2.remove or g2.rescue) and g2.alternative_transcript_rescue != set():
                            if bool(g.alternative_transcript_rescue.intersection(g2.alternative_transcript_rescue)):
                                temp_group = g.alternative_transcript_rescue | g2.alternative_transcript_rescue
                                temp_group.add(g.id)
                                temp_group.add(g2.id)
                                found = False
                                index = 0
                                for n, group in enumerate(gene_groups):
                                    if bool(temp_group.intersection(group)):
                                        found = True
                                        index = n
                                        break
                                if found:
                                    gene_groups[index] = gene_groups[index] | temp_group
                                else:
                                    gene_groups.append(temp_group)

        merge_genes = set()
        for gene_set in gene_groups:
            merge_genes = merge_genes | gene_set
            best_reliable = ""
            best_unreliable = ""
            for g_id in gene_set:
                chrom = self.all_gene_ids[g_id]
                if not self.chrs[chrom][g_id].remove or self.chrs[chrom][g_id].rescue:
                    best_reliable = g_id
                else:
                    best_unreliable = g_id
            
            for g_id in gene_set:
                chrom = self.all_gene_ids[g_id]
                if not self.chrs[chrom][g_id].remove or self.chrs[chrom][g_id].rescue:
                    query_best = self.chrs[chrom][g_id].compare_protein_blast_hits(self.chrs[chrom][best_reliable], source_priority)
                    if query_best:
                        best_reliable = g_id
                else:
                    query_best = self.chrs[chrom][g_id].compare_protein_blast_hits(self.chrs[chrom][best_unreliable], source_priority)
                    if query_best:
                        best_unreliable = g_id

            for g_id in gene_set:
                chrom = self.all_gene_ids[g_id]

                if g_id == best_reliable:
                    continue
                elif g_id == best_unreliable:
                    for t in self.chrs[chrom][best_unreliable].transcripts.values():
                        if t.main:
                            t_copy = t.copy()
                            t_copy.id = "alternative_transcript"
                            t_copy.parents = [g_id]
                            t_copy.symbols = []
                            t_copy.names = []
                            t_copy.synonyms = []
                    self.chrs[chrom][best_reliable].transcripts["alternative_transcript"] = t_copy.copy()
                    del t_copy
                else:
                    self.chrs[chrom][g_id].rescue = False
                    self.chrs[chrom][g_id].remove = True

        for genes in self.chrs.values():
            for g in genes.values():
                g.alternative_transcript_rescue = list(g.alternative_transcript_rescue)
                progress_bar.update(1)
                if g.id not in merge_genes and g.alternative_transcript_rescue != []:
                    best = g.alternative_transcript_rescue[0]
                    for alt in g.alternative_transcript_rescue:
                        if alt == best:
                            continue
                        query_best = g.compare_protein_blast_hits(self.chrs[g.ch][alt], source_priority)
                        if not query_best:
                            best = alt

                    for t in self.chrs[g.ch][best].transcripts.values():
                        if t.main:
                            t_copy = t.copy()
                            t_copy.id = "alternative_transcript"
                            t_copy.parents = [g.id]
                            t_copy.symbols = []
                            t_copy.names = []
                            t_copy.synonyms = []
                    g.transcripts[t_copy.id] = t_copy.copy()
                    del t_copy
        progress_bar.close()
        self.update()

    def remove_exon_overlaps(self, source_priority, blast:bool=False):
        for genes in self.chrs.values():
            for g in genes.values():
                if not g.remove or g.rescue:
                    for o in g.overlaps["self"]:
                        if (not self.chrs[g.ch][o.id].remove or self.chrs[g.ch][o.id].rescue) and (not g.remove or g.rescue):
                            if o.exon_query_percent > 0 or o.exon_target_percent > 0:
                                if blast:
                                    query_best = g.compare_protein_blast_hits(self.chrs[g.ch][o.id], source_priority)
                                else:
                                    query_best = g.longer_CDS(self.chrs[g.ch][o.id])
                                if query_best:
                                    self.chrs[g.ch][o.id].remove = True
                                    self.chrs[g.ch][o.id].rescue = False
                                else:
                                    g.remove = True
                                    g.rescue = False

    def remove_UTRs_from_exon_overlaps(self):

        for genes in self.chrs.values():
            for g in genes.values():
                if not g.remove or g.rescue:
                    for o in g.overlaps["self"]:
                        if (not self.chrs[g.ch][o.id].remove or self.chrs[g.ch][o.id].rescue) and (not g.remove or g.rescue):
                            if o.exon_query_percent > 0 or o.exon_target_percent > 0:
                                self.chrs[g.ch][o.id].clear_UTRs()
                                g.clear_UTRs()

        self.sorted = False

    def remove_CDS_overlaps(self, source_priority, blast:bool=False, anti:bool=True):
        for genes in self.chrs.values():
            for g in genes.values():
                if not g.remove or g.rescue:
                    for o in g.overlaps["self"]:
                        overlapping_CDS = False
                        if anti:
                            if o.score >= 5 or (o.score == 1 and o.antiscore >= 5):
                                overlapping_CDS = True
                        else:
                            if o.score >= 5:
                                overlapping_CDS = True

                        if overlapping_CDS:
                            if (not self.chrs[g.ch][o.id].remove or self.chrs[g.ch][o.id].rescue) and (not g.remove or g.rescue):
                                if blast:
                                    query_best = g.compare_protein_blast_hits(self.chrs[g.ch][o.id], source_priority)
                                else:
                                    query_best = g.longer_CDS(self.chrs[g.ch][o.id])
                                if query_best:
                                    self.chrs[g.ch][o.id].remove = True
                                    self.chrs[g.ch][o.id].rescue = False
                                else:
                                    g.remove = True
                                    g.rescue = False

    def remove_fully_intron_nested_genes(self):
        for genes in self.chrs.values():
            for g in genes.values():
                if not g.remove or g.rescue:
                    if g.intron_nested and not g.UTR_intron_nested:
                        for t in g.transcripts.values():
                            if t.main:
                                for c in t.CDSs.values():
                                    if c.main:
                                        if c.size <= 450 or t.coding_ratio < 0.4:
                                            g.remove = True
                                            g.rescue = False

    def rescue_longer_same_frame_CDS(self, reliable_sources:list=["AUGUSTUS", "GeneMark.hmm3"]):

        for genes in self.chrs.values():
            for g in genes.values():
                if not g.remove or g.rescue:

                    for t in g.transcripts.values():
                        if t.main:
                            for c in t.CDSs.values():
                                if c.main:
                                    main_CDS_size = c.size

                    posible_alternative_transcripts = []

                    for o in g.overlaps["self"]:
                        if o.CDSs_in_both:

                            for t in self.chrs[g.ch][o.id].transcripts.values():
                                if t.main:
                                    for c in t.CDSs.values():
                                        if c.main:
                                            overlap_main_CDS_size = c.size

                            if (self.chrs[g.ch][o.id].source in reliable_sources) and ((o.full_protein_overlaps >= 2) or (o.protein_query_percent >= 90)) and (overlap_main_CDS_size > main_CDS_size):

                                for t in self.chrs[g.ch][o.id].transcripts.values():
                                    if t.main:
                                        t_copy = t.copy()
                                        t_copy.id = "alternative_transcript2"
                                        t_copy.parents = [g.id]
                                        t_copy.symbols = []
                                        t_copy.names = []
                                        t_copy.synonyms = []

                                posible_alternative_transcripts.append(t_copy.copy())
                                del t_copy

                    best_candidate_size = 0
                    best_candidate = None

                    if posible_alternative_transcripts:

                        for t_candidate in posible_alternative_transcripts:

                            if t_candidate.size > best_candidate_size:

                                best_candidate = t_candidate
                                best_candidate_size = t_candidate.size

                        g.transcripts[best_candidate.id] = best_candidate.copy()

        self.update()

    def make_alternative_genes_into_transcripts(self):

        correspondence = {}

        to_remove = []
        for chrom, genes in self.chrs.items():
            for g in genes.values():
                if (not g.remove or g.rescue) and (g.id not in correspondence):
                    for o in g.overlaps["self"]:
                        if o.score >= 5:

                            if o.id in correspondence:

                                correspondence_gene = correspondence[o.id]
                                correspondence[g.id] = correspondence_gene
                                to_remove.append((g.id, chrom))
                                for t in self.chrs[chrom][g.id].transcripts.values():
                                    self.chrs[chrom][correspondence_gene].transcripts[t.id] = t.copy()

                            else:

                                correspondence[o.id] = g.id
                                to_remove.append((o.id, chrom))
                                for t in self.chrs[chrom][o.id].transcripts.values():
                                    g.transcripts[t.id] = t.copy()

        for g, chrom in to_remove:
            if g in self.chrs[chrom]:
                del self.chrs[chrom][g]

        self.update()

    def find_best_gene_model_exon_num_overlaps(self, source_priority, blast:bool=False, exon_num:int=2):
        """
        For genes with more than X exons exactly the same.
        """
        for genes in self.chrs.values():
            for g in genes.values():
                if not g.remove or g.rescue:
                    for o in g.overlaps["self"]:
                        if (not self.chrs[g.ch][o.id].remove or self.chrs[g.ch][o.id].rescue) and (not g.remove or g.rescue):
                            if o.full_exon_overlaps >= exon_num:
                                if blast:
                                    query_best = g.compare_protein_blast_hits(self.chrs[g.ch][o.id], source_priority)
                                else:
                                    query_best = g.longer_CDS(self.chrs[g.ch][o.id])
                                if query_best:
                                    self.chrs[g.ch][o.id].remove = True
                                    self.chrs[g.ch][o.id].rescue = False
                                else:
                                    g.remove = True
                                    g.rescue = False

    def find_best_gene_model_nested_overlaps(self, source_priority, blast=False):
        """
        For genes fully contained in other genes which have exon overlap choose best blast hit.
        """
        for genes in self.chrs.values():
            for g in genes.values():
                if not g.remove or g.rescue:
                    for o in g.overlaps["self"]:
                        if (not self.chrs[g.ch][o.id].remove or self.chrs[g.ch][o.id].rescue) and (not g.remove or g.rescue):
                            if o.gene_query_percent >= 100 or o.gene_target_percent >= 100:
                                if o.exon_query_percent > 0:
                                    if blast:
                                        query_best = g.compare_protein_blast_hits(self.chrs[g.ch][o.id], source_priority)
                                    else:
                                        query_best = g.longer_CDS(self.chrs[g.ch][o.id])
                                    if query_best:
                                        self.chrs[g.ch][o.id].remove = True
                                        self.chrs[g.ch][o.id].rescue = False
                                    else:
                                        g.remove = True
                                        g.rescue = False

    def gene_count(self):
        gene_objects = 0
        unique_gene_ids_in_overlaps = set()
        for genes in self.chrs.values():
            gene_objects += len(genes)
            for g in genes.values():
                for o in g.overlaps["self"]:
                    unique_gene_ids_in_overlaps.add(o.id)
        print(f"There are {gene_objects} gene objects and {len(self.all_gene_ids)} genes in all gene ids and {len(unique_gene_ids_in_overlaps)} ids contained in overlaps.")

    def remove_redundancy(self, source_priority:list, hard_masked_genome:object):
        self.remove_duplicate_transcripts()
        self.make_alternative_transcripts_into_genes()
        self.detect_gene_overlaps()
        self.calculate_transcript_masking(hard_masked_genome=hard_masked_genome)
        self.mark_noisy_genes()
        self.remove_genes()
        self.mark_transcriptomic_supported_genes()
        self.mark_abinitio_supported_genes()
        self.add_reliable_CDS_evidence_score()
        self.find_best_gene_model(source_priority)
        self.mark_overlap_with_reliable_genes()
        self.find_best_gene_model(source_priority, just_with_reliables=False)

        self.add_better_ab_initio_models_as_alternative_transcripts(source_priority, reliable_sources=["AUGUSTUS", "Liftoff", "GeneMark.hmm3"])
        self.detect_gene_overlaps()

        self.rescue_longer_same_frame_CDS()
        self.detect_gene_overlaps()

        self.remove_CDS_overlaps(source_priority)
        self.mark_intron_nesting()
        self.remove_fully_intron_nested_genes()

        self.make_alternative_transcripts_into_genes()
        self.detect_gene_overlaps()

        self.mark_overlap_with_other_selected_CDSs()
        self.mark_overlap_with_other_selected_exons()
        self.select_best_possible_non_overlapping_UTR()
        self.detect_gene_overlaps()

        self.mark_overlap_with_other_selected_CDSs()
        self.mark_overlap_with_other_selected_exons()
        self.select_best_possible_non_overlapping_UTR(exon=True)

        self.remove_genes()
        self.detect_gene_overlaps()

        self.make_alternative_genes_into_transcripts()
        self.detect_gene_overlaps()
        self.find_best_gene_model_nested_overlaps(source_priority)
        self.find_best_gene_model_exon_num_overlaps(source_priority)
        self.remove_exon_overlaps(source_priority)
        self.remove_UTRs_from_exon_overlaps()
        self.remove_genes()
        self.update()
        self.detect_gene_overlaps()

    def remove_genes_with_small_CDSs(self, CDS_threshold:int=200):

        removed_any = False

        for genes in self.chrs.values():
            for g in genes.values():
                remove = False
                for t in g.transcripts.values():
                    if t.main:
                        for c in t.CDSs.values():
                            if c.main:
                                if c.size < CDS_threshold:
                                    remove = True
                                    removed_any = True
                if remove:
                    g.rescue = False
                    g.remove = True
        self.remove_genes()

        if removed_any:
            self.small_cds_removed = True
        self.update(rename_ids_minimal=True)

    def remove_TE_genes(self):

        removed_any = False

        for genes in self.chrs.values():
            for g in genes.values():
                if g.transposable:
                    g.rescue = False
                    g.remove = True
                    removed_any = True

        self.remove_genes()

        if removed_any:
            self.transposable_removed = True
        self.update(rename_ids_minimal=True)

    def remove_non_TE_genes(self):

        removed_any = False

        for genes in self.chrs.values():
            for g in genes.values():
                if not g.transposable:
                    g.rescue = False
                    g.remove = True
                    removed_any = True

        self.remove_genes()
        
        if removed_any:
            self.non_transposable_removed = True

        self.update(rename_ids_minimal=True)

    def remove_non_coding_genes_and_transcripts(self):

        removed_any = False

        for genes in self.chrs.values():
            for g in genes.values():
                if g.noncoding == True and g.coding == False:
                    g.rescue = False
                    g.remove = True
                    removed_any = True

        self.remove_genes()

        removed_any = self.remove_non_coding_transcripts_from_coding_genes(removed_any, False)

        if removed_any:
            self.non_coding_removed = True

        self.update(rename_ids_minimal=True)

    def remove_coding_genes_and_transcripts(self):

        removed_any = False

        for genes in self.chrs.values():
            for g in genes.values():
                if g.noncoding == False and g.coding == True:
                    g.rescue = False
                    g.remove = True
                    removed_any = True

        self.remove_genes()

        removed_any = self.remove_coding_transcripts_from_non_coding_genes(removed_any, False)

        if removed_any:
            self.coding_removed = True

        self.update(rename_ids_minimal=True)

    def remove_coding_transcripts_from_non_coding_genes(self, removed_any:bool=False, update=True):
        transcripts_to_remove = []

        for chrom, genes in self.chrs.items():
            for g in genes.values():
                if g.coding:
                    for t in g.transcripts.values():
                        if t.coding == True:
                            transcripts_to_remove.append((chrom, g.id, t.id))
        
        for chrom, g_id, t_id in transcripts_to_remove:
            del self.chrs[chrom][g_id].transcripts[t_id]

        if transcripts_to_remove != []:
            removed_any = True

        if removed_any:
            self.coding_removed = True
        
        self.clear_overlaps()

        if update:
            self.update(rename_ids_minimal=True)

        return removed_any

    def remove_non_coding_transcripts_from_coding_genes(self, removed_any:bool=False, update=True):
        transcripts_to_remove = []

        for chrom, genes in self.chrs.items():
            for g in genes.values():
                if g.noncoding:
                    for t in g.transcripts.values():
                        if t.coding == False:
                            transcripts_to_remove.append((chrom, g.id, t.id))
        
        for chrom, g_id, t_id in transcripts_to_remove:
            del self.chrs[chrom][g_id].transcripts[t_id]

        if transcripts_to_remove != []:
            removed_any = True

        if removed_any:
            self.non_coding_removed = True

        self.clear_overlaps()

        if update:
            self.update(rename_ids_minimal=True)

        return removed_any

    def clear_gene_names_and_symbols(self):
        for genes in self.chrs.values():
            for g in genes.values():
                g.names = []
                g.symbols = []
                g.synonyms = []
        self.update()

    def remove_genes_without_symbols(self):
        for genes in self.chrs.values():
            for g in genes.values():
                if g.symbols == []:
                    g.remove = True
        self.remove_genes()
        self.update()

    def rename_chromosomes(self, equivalences, dap:bool=False):

        renamed_scaffolds = False
        for old, new in equivalences.items():
            if new != old and old in self.chrs:
                renamed_scaffolds = True

        self.chrs = {equivalences.get(k, k): v for k, v in self.chrs.items()}

        for chrom, genes in self.chrs.items():
            for g in genes.values():
                g.ch = chrom
                for t in g.transcripts.values():
                    t.ch = chrom
                    if hasattr(t, "promoter"):
                        t.promoter.ch = chrom
                    for c in t.CDSs.values():
                        c.ch = chrom
                        for cs in c.CDS_segments:
                            cs.ch = chrom
                        if hasattr(c, "protein"):
                            if c.protein != None:
                                c.protein.ch = chrom
                        for u in c.UTRs:
                            u.ch = chrom
                    for e in t.exons:
                        e.ch = chrom
                    for i in t.introns:
                        i.ch = chrom

        for o in self.atypical_features:
            if o.ch in equivalences:
                o.ch = equivalences[o.ch]

        if dap:
            self.dapfit = True
            if renamed_scaffolds:
                self.dapmod = True
        elif renamed_scaffolds:
            self.confrenamed = True

        self.update_suffixes()

    def add_gene_symbols(self, file_path:str, just_gene_names:bool=True, clear:bool=True):
        if clear:
            self.clear_gene_names_and_symbols()

        if file_path.endswith(".xlsx"):
            df = pd.read_excel(file_path, skiprows=1, dtype=str)
        else:
            df = pd.read_csv(file_path, sep="\t", dtype=str)
        df = df.fillna("")

        pseudogene_col_exists = False
        synonym_col_exists = False

        if "pseudogene" in df.columns:
            pseudogene_col_exists = True
        
        if "gene name synonym(s)" in df.columns:
            synonym_col_exists = True

        for index, row in df.iterrows():
            gene_name = df.iloc[index, 0]
            gene_id = df.iloc[index, 1]
            if pseudogene_col_exists:
                pseudogene = row["pseudogene"]
            if synonym_col_exists:
                gene_synonyms = row["gene name synonym(s)"].split("; ")
            ch = self.all_gene_ids[gene_id]
            self.chrs[ch][gene_id].symbols.append(gene_name)
            if not just_gene_names:
                if synonym_col_exists:
                    self.chrs[ch][gene_id].synonyms += gene_synonyms
                if pseudogene_col_exists:
                    if pseudogene == "pseudogene" or pseudogene == "Pseudogene":
                        self.chrs[ch][gene_id].pseudogene = True
                    elif pseudogene == "gene" or pseudogene == "Gene":
                        self.chrs[ch][gene_id].pseudogene = False

        self.symbols_added = True

        self.update_attributes(clean=True, aliases=False, extra_attributes=False, symbols=True, symbols_as_descriptors=True)

    def release(self, name, id, source_name, id_prefix, spacer:int=10, suffix:str="", custom_path:str="", tag:str=".gff3", skip_atypical_fts:bool=True, main_only:bool=False, UTRs:bool=True, clear_aliases=True, extra_attributes=False):
        if clear_aliases:
            self.clear_aliases()
        self.name = name
        self.id = id
        for genes in self.chrs.values():
            for g in genes.values():
                g.source = source_name
                for t in g.transcripts.values():
                    t.source = source_name
                    for e in t.exons:
                        e.source = source_name
                    for c in t.CDSs.values():
                        c.source = source_name
                        for cs in c.CDS_segments:
                            cs.source = source_name
                        for u in c.UTRs:
                            u.source = source_name

        self.update(extra_attributes=extra_attributes)
        self.rename_ids(prefix=id_prefix, spacer=spacer, suffix=suffix)
        self.update(extra_attributes=extra_attributes)
        self.export_gff(custom_path=custom_path, tag=tag, skip_atypical_fts=skip_atypical_fts, main_only=main_only, UTRs=UTRs)

    def __str__(self):
        return str(self.id)