import time
import textwrap
import pandas as pd
from os import system
from Bio import SeqIO
from pathlib import Path

class Genome():
    def __init__(self, name:str, genome_file_path:str, chromosome_dict:dict={}):
        start = time.time()
        self.name = name
        self.file = genome_file_path
        self.path = "/".join(self.file.split("/")[0:-1]) + "/"
        # Creating a dictionary with the genome sequence of each chromosome or scaffold, still referred as chromosomes in the code
        self.chseqs = {}

        with open (self.file, "r", encoding="utf-8") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                ch = record.id
                if ch in self.chseqs:
                    print((f"Error: chromosome/scaffold {ch} is repeated in {self.name}, genome (file: {self.file})"))
                self.chseqs[ch] = str(record.seq).upper()

        self.features = list(self.chseqs.keys())
        self.size = sum(len(seq) for seq in self.chseqs.values())

        self.dapfit = True
        for chrom in self.chseqs:
            if not chrom.startswith("chr") and not chrom.startswith("Chr"):
                self.dapfit = False
            else:
                try:
                    number = int(chrom.split("r")[1])
                except:
                    self.dapfit = False

        self.scaffolds = False
        for chrom in self.chseqs:
            if chrom.lower() == "chr00" or chrom.lower() == "ch00" or chrom.lower() == "chrun" or chrom.lower() == "chun":
                self.scaffold = True
            if chrom in chromosome_dict:
                continue
            if not chrom.startswith("ch") and not chrom.startswith("Ch"):
                self.scaffolds = True      

        self.dapmod = False

        if "_dapmod" in genome_file_path:
            self.dapmod = True

        now = time.time()
        lapse = now - start
        print(f"\n{self.name} genome chromosomes/scaffolds: {self.features[:30]} ...")
        print(f"\nCreating {self.name} genome object took {round(lapse, 1)} seconds\n")

    def export_feature_sizes(self, custom_path:str=""):

        if custom_path == "":
            export_folder = self.path + "genome_feature_sizes/"
        else:
            export_folder = custom_path
        if export_folder[-1] != "/":
            export_folder += "/"
        system(f"mkdir -p {export_folder}")

        if not self.dapmod:
            tag = f"{self.name}_genome_feature_sizes.tsv"
        else:
            tag = f"{self.name}_genome_feature_sizes_dapmod.tsv"

        out = ""
        for chrom, seq in self.chseqs.items():
            out += f"{chrom}\t{len(seq)}\n"
        f_out = open(f"{export_folder}{tag}", "w", encoding="utf-8")
        f_out.write(out)
        f_out.close()
    
    def rename_features_dap(self, output_folder:str="", return_equivalences:bool=False, export:bool=False, chromosome_dict:dict={}):
        """
        Renames scaffolds and chromosomes for DAP-Seq.
        """

        if not self.dapfit:

            equivalences = {}

            scaffold_count = 0
            for feature in self.chseqs:
                if feature in chromosome_dict:
                    equivalences[feature] = chromosome_dict[feature]
                elif feature.startswith("chr") or feature.startswith("Chr"):
                    try:
                        number = "{:02d}".format(int(feature.split("r")[1]))
                        equivalences[feature] = f"chr{number}"
                    except:
                        scaffold_count += 1
                        number = "{:08d}".format(scaffold_count)
                        equivalences[feature] = f"chr{number}"
                elif feature.startswith("ch") or feature.startswith("Ch"):
                    try:
                        number = "{:02d}".format(int(feature.split("h")[1]))
                        equivalences[feature] = f"chr{number}"
                    except:
                        scaffold_count += 1
                        number = "{:08d}".format(scaffold_count)
                        equivalences[feature] = f"chr{number}"
                else:
                    scaffold_count += 1
                    number = "{:08d}".format(scaffold_count)
                    equivalences[feature] = f"chr{number}"

            self.chseqs = {equivalences.get(k, k): v for k, v in self.chseqs.items()}
            self.dapmod = True
            self.dapfit = True

            if export:
                self.export(output_folder=output_folder, file=".fasta")

            if return_equivalences:
                return equivalences
            
        else:

            print(f"Warning: {self.name} genome is already fit for DAP-Seq analysis, so it won't be modified.")

    def remove_scaffolds(self, output_folder:str="", export:bool=False, chromosome_dict:dict={}):
        if self.scaffolds:
            new_chseqs = {}

            for scaffold, values in self.chseqs.items():
                if scaffold.lower() == "chr00" or scaffold.lower() == "ch00" or scaffold.lower() == "chrun" or scaffold.lower() == "chun":
                    continue
                if scaffold in chromosome_dict or scaffold.startswith("Ch") or scaffold.startswith("ch"):
                    new_chseqs[scaffold] = values

            self.chseqs = new_chseqs.copy()
            del new_chseqs
            self.scaffolds = False

            if export:
                self.export(output_folder=output_folder, file=".fasta")

    def export(self, output_folder:str="", file:str=".fasta"):
        if not output_folder:
            export_folder = self.path + "out_genomes/"
        else:
            export_folder = output_folder
        if export_folder[-1] != "/":
            export_folder += "/"

        system(f"mkdir -p {export_folder}")

        if file == ".fasta":
            added = ""
            if self.dapmod:
                added += "_dapmod"
            elif self.dapfit:
                added += "_dapfit"
            if not self.scaffolds:
                added += "_scffree"

        file = f"{self.name}{added}{file}"

        out = ""

        for feature, seq in self.chseqs.items():
            out += f">{feature}\n{seq}\n"

        f_out = open(f"{export_folder}{file}", "w", encoding="utf-8")
        f_out.write(out)
        f_out.close()
        
    def extract_peak_sequences(self, output_file_name, DAPseq_output_file, output_folder:str="", top=600):

        export_folder = Path(output_folder or self.path + "out_peak_seqs/")
        export_folder.mkdir(parents=True, exist_ok=True)

        df = pd.read_csv(DAPseq_output_file, delimiter='\t', dtype=str)
        df.dropna(how='all', inplace=True)

        if top != 'all':
            top_df = df
        else:
            top_df = df.nlargest(top, 'peak signal value')

        peaks = {
            f"{row['seqnames']}_{row['feature']}_{int(row['peak'].split('.')[1])-100}:{int(row['peak'].split('.')[1])+100}":
            self.chseqs[row['seqnames']][int(row['peak'].split('.')[1])-100:int(row['peak'].split('.')[1])+100]
            for _, row in top_df.iterrows()
        }

        with open(export_folder / output_file_name, "w", encoding="utf-8") as f_out:
            for header, seq in peaks.items():
                f_out.write(f'>{header}\n')
                f_out.write(f'{textwrap.fill(seq, width=60)}\n')
