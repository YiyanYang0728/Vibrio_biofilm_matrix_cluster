gene=$1
gid=$(echo ${gene} | cut -f1 -d.)
phage_annot=phage_pharokka/${gid}*/*.tsv
ctg_start=`dirname ${phage_annot} | awk -F"_" '{split($NF, a, "-");print a[1]}'`
gene_gff=/gpfs/gsfs12/users/yangy34/projects/vibrio_YJ_v5/VPS_cluster/genomes/gff/${gid}.*.gff
# echo $ctg_start

grep -w ${gene} ${gene_gff} | cut -f1,4,5,7|awk -F"\t" -v ID=${gene} 'BEGIN{OFS="\t"}{print $1,$2,$3,ID,"-",$4}' |bedtools sort > ${gene}_gff.bed
ctg=`cut -f1 ${gene}_gff.bed|sort -u`
cut -f1,2,3,4 ${phage_annot} | awk -v ctg=${ctg} -v ctg_start=${ctg_start} 'BEGIN{OFS="\t"}NR>1{ID=$1;$1=ctg;if($2>$3){print $1,$3+ctg_start-1,$2+ctg_start-1,ID,"-",$4}else{print $1,$2+ctg_start-1,$3+ctg_start-1,ID,"-",$4}}' | bedtools sort > ${gene}_prophage.bed

bedtools closest -a ${gene}_gff.bed -b ${gene}_prophage.bed -s -d | awk -F"\t" -v ctg=${ctg} '$NF==0{print $4"\t"$10"\t"ctg}'
rm ${gene}_gff.bed ${gene}_prophage.bed