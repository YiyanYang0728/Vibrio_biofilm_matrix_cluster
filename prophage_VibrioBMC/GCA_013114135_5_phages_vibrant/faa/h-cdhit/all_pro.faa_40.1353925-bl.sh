#!/bin/sh
#PBS -v PATH
#$ -v PATH


para=$1
cd /gpfs/gsfs12/users/yangy34/projects/vibrio_paper/prophage_VibrioBMC/GCA_013114135_5_phages_vibrant/faa/h-cdhit
./all_pro.faa_40.1353925-bl.pl 0 $para &
wait

