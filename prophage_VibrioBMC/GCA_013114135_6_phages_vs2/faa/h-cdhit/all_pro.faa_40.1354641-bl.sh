#!/bin/sh
#PBS -v PATH
#$ -v PATH


para=$1
cd /gpfs/gsfs12/users/yangy34/projects/vibrio_paper/prophage_VibrioBMC/GCA_013114135_6_phages_vs2/faa/h-cdhit
./all_pro.faa_40.1354641-bl.pl 0 $para &
wait

