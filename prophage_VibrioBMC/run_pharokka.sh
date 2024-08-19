#!/bin/bash
#SBATCH --partition norm
#SBATCH --job-name pharokka
#SBATCH --cpus-per-task=32
#SBATCH --mem=50g
#SBATCH --time=5-00:00:00
source /data/$USER/conda/etc/profile.d/conda.sh
conda activate /gpfs/gsfs12/users/yangy34/softwares/mambaforge/envs/pharokka-env
infile=$1
id=$(grep ">" ${infile}|head -n1|sed "s/>//g")
pharokka.py -f -i ${infile} -o phage_pharokka/${id} -d /data/yangy34/database/pharokka -t 32
