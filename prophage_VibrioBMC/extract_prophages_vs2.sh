gid=$1
outdir=${2:-$PWD}
# get prophage numbers & get contig, start and end pos for each phage
file=/data/yangy34/projects/Vibrio_phages/run_vs2/vs2_out/${gid}*/final-viral-boundary.tsv
tail -n+2 ${file} | cut -f1,4,5,17,18,29,31 > ${outdir}/${gid}_prophages_vs2.tsv
ct=`cat ${outdir}/${gid}_prophages_vs2.tsv|wc -l`
echo -e "$gid\t$ct"
#rm ${outdir}/${gid}_prophages_vs2.tsv
