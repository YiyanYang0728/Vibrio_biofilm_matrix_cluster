#!/bin/bash
#SBATCH --partition quick
#SBATCH --job-name mafft
#SBATCH --cpus-per-task=64
#SBATCH --mem=50g
#SBATCH --time=04:00:00
source /data/$USER/conda/etc/profile.d/conda.sh
conda activate /data/yangy34/conda/envs/bio-env
cd /gpfs/gsfs12/users/yangy34/projects/vibrio_paper/Step3.annotate_RbmB_RbmA
#mafft --dash --thread 32 --localpair ${f} &> out/${prefix}.out
mafft --thread 64 --maxiterate 3 --localpair rhbh_domain_len_filtered_rmdup.fasta > rhbh_domain_len_filtered_rmdup.mafft.aln.4
