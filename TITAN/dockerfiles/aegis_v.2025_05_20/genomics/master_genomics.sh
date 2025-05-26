#!/bin/bash
#SBATCH --job-name=genomics       # Job name (showed with squeue)
#SBATCH --output=genomics_output_%j.log  # Standard output and error log
#SBATCH --qos=medium                 # QoS: short,medium,long,long-mem
#SBATCH --nodes=1                   # Required only 1 node
#SBATCH --ntasks=1                  # Required only 1 task
#SBATCH --cpus-per-task=2           # Required only 1 cpu
#SBATCH --mem=400G                   # Required 10GB of memory
#SBATCH --time=2-00:00:00             # Required 5 minutes of execution time

module load python/3.11.4
cd /storage/tom/repositorios_antonio/genomics
python3 -u test_new_annotationClass_improved_garnatxa.py
python3 /storage/tom/Antonio_TFM/telegram_bot.py