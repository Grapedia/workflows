#!/bin/bash

#SBATCH --job-name=genomics
#SBATCH -o logs/annotation_objects/annotation_objects_manager_%j_out.log
#SBATCH -e logs/annotation_objects/annotation_objects_manager_%j_err.log
#SBATCH --ntasks=1    
#SBATCH --time=3-00:00:0
#SBATCH --mem=110gb
#SBATCH --cpus-per-task=2
#SBATCH --qos=long

module load python/3.11.4


python3 -u annotation_objects.py --combine_transcripts True --promoters True --species cannabis

#python3 -u annotation_objects.py --combine_transcripts True --promoters True --species peach
#example uses:

# will generate all annotation objects for all species

#python3 -u annotation_objects.py

#python3 -u annotation_objects.py --species cannabis --annotations CBDRx_on_CBDRx PurpleKush_on_PurpleKush

#python3 -u annotation_objects.py --combine_transcripts True --promoters True --species PN40024_T2T chickpea tomato hop snapdragon chlamydomonas mulberry