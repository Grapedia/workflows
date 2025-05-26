# genomes-annotations


For now just some classes defined to import genomes from fasta files and annotations from gff3 files. Some general bioinformatics functions are also included as a separate module. Use the overlap function to find correspondences between different annotation files associated to the same genome. If one of those files comes from a liftoff and is indicated when loading the annotation file, synteny conservation during liftoff is also worked out and present in the overlap function output.

requirements
python 3.11.4
kaleido                       0.2.1
pandas                        1.1.5
plotly                        5.18.0
biopython                     1.79