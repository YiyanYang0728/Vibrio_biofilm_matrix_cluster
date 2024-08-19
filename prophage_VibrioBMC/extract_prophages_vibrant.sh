gid=$1
outdir=${2:-$PWD}

doc=`echo ${gid}|awk '{print substr($0,1,3)"/"substr($0,5,3)"/"substr($0,8,3)"/"substr($0,11,3)}'`
gid_1=`echo ${gid} | cut -f1 -d.`
tail -n+2 /data/yangy34/projects/Vibrio_phages/data/${doc}/${gid_1}_vibrant/VIBRANT_results_${gid_1}/VIBRANT_integrated_prophage_coordinates_${gid_1}.tsv|cut -f1,6,7,8 > ${outdir}/${gid}_prophages_vibrant.tsv
# get prophage numbers & get contig, start and end pos for each phage
ct=`cat ${outdir}/${gid}_prophages_vibrant.tsv|wc -l`
echo -e "$gid\t$ct"
#rm ${outdir}/${gid}_prophages_vibrant.tsv
