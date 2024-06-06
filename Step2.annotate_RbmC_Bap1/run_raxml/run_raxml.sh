#!/bin/bash
#SBATCH --partition norm
#SBATCH --job-name raxml
#SBATCH --cpus-per-task=56
#SBATCH --mem=50g
#SBATCH --time=5-00:00:00
source /data/$USER/conda/etc/profile.d/conda.sh
conda activate /data/yangy34/conda/envs/bio-env

cd /data/yangy34/projects/vibrio_paper/process/Step2.annotate_RbmC_Bap1/run_raxml
raxmlHPC-PTHREADS-SSE3 -T 56 -m GTRGAMMA -p 12345 -q dna12_3.partition.txt -s all_genes.rmdup.codon.aln.fa -n raxml.out