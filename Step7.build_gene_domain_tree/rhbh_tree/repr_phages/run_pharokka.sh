#!/bin/bash
#SBATCH --partition norm
#SBATCH --job-name pharokka
#SBATCH --cpus-per-task=32
#SBATCH --mem=50g
#SBATCH --time=5-00:00:00
source /data/$USER/conda/etc/profile.d/conda.sh
conda activate /gpfs/gsfs12/users/yangy34/softwares/mambaforge/envs/pharokka-env
cd /gpfs/gsfs12/users/yangy34/projects/vibrio_paper/Step7.build_gene_domain_tree/rhbh_tree/repr_phages
infile=$1
id=`grep ">" ${infile}|sed "s/>//g;s/-/_/g"`
#echo $id
pharokka.py -f  -i ${infile} -o phage_pharokka/${id} -d /data/yangy34/database/pharokka -t 32
