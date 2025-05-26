import time
import os
from modules.tools import pickle_load, pickle_save
from modules.annotation import Annotation
from modules.genome import Genome


gff_file = "../../genomes_and_annotation/genemark.gff3"
pickle_tag = "../../genomes_and_annotation/genemark"
riesling_genome = "../../genomes_and_annotation/riesling.hap1.chromosomes.phased.fa"
location = f"{pickle_tag}_annotation.pkl"

genome = Genome("genemark", riesling_genome)

if os.path.isfile(location):
    a = pickle_load(location)
else:
    a = Annotation("genemark", gff_file, genome)
    pickle_save(location, a)

a.update_stats(custom_path="../../genomes_and_annotation/", export=True, genome=genome)