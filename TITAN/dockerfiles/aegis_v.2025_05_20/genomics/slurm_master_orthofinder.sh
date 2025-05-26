#!/bin/bash

#SBATCH --job-name=orthofinder
#SBATCH -o logs/liftoffs/orthofinder_%j_out.log
#SBATCH -e logs/liftoffs/orthofinder_%j_err.log
#SBATCH --ntasks=1    
#SBATCH --time=20:00:00       
#SBATCH --mem=20gb
#SBATCH --cpus-per-task=8
#SBATCH --qos=short

module load python/3.11.4

python3 /storage/tom/apps/OrthoFinder_source/orthofinder.py -f /storage/tom/genomes_and_annotation/cannabis/CBDRx_annotation_transfer/orthofinder_jamaicanlion

python3 /storage/tom/apps/OrthoFinder_source/orthofinder.py -f /storage/tom/genomes_and_annotation/cannabis/CBDRx_annotation_transfer/orthofinder_finola

python3 /storage/tom/apps/OrthoFinder_source/orthofinder.py -f /storage/tom/genomes_and_annotation/cannabis/CBDRx_annotation_transfer/orthofinder_purplekush

