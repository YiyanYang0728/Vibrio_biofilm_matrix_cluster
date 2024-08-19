gene=$1
gid=`echo ${gene}|cut -f1 -d. `
/gpfs/gsfs12/users/yangy34/projects/vibrio_YJ_v5/VPS_cluster/genomes/RbmB_genes/Vibrio_rhbh_evol/get_gene_context.sh ${gene} 0 > ${gene}.tmp
ctg=`cut -f1 ${gene}.tmp`
start=`cut -f2 ${gene}.tmp`
end=`cut -f3 ${gene}.tmp`

# doc=`echo ${gid}|awk '{print substr($0,1,3)"/"substr($0,5,3)"/"substr($0,8,3)"/"substr($0,11,3)}'`
fgrep -w ${ctg} /data/yangy34/projects/vibrio_paper/Step8.find_prophage_regions/vs2_prophage_regions.tsv > ${gene}.vs2.tmp
printout=0
if [ -s ${gene}.vs2.tmp ]
then

    while IFS=$'\t' read -r p_ctg ns ne; do
        if [ "${p_ctg}" = "${ctg}" ] && [ "${ns}" -le "${start}" ] && [ "${ne}" -ge "${end}" ]; then
            echo -e "${gene}\tIn_prophage\tThis gene ${ctg}:${start}-${end} is in prophage ${ctg}:${ns}-${ne}"
            printout=1
            break
        fi
    done < ${gene}.vs2.tmp

    if [ $printout -eq 0 ]; then
        while IFS=$'\t' read -r p_ctg ns ne; do
            if [ "${p_ctg}" = "${ctg}" ]; then
                echo -e "${gene}\tSame_contig\tThis gene ${ctg}:${start}-${end} is close to the prophage ${ctg}:${ns}-${ne}"
                printout=1
                break
            fi
        done < ${gene}.vs2.tmp
    fi

else
    echo -e "${gene}\tNot_in_prophage\tNo prophage in this contig"
fi
rm ${gene}.vs2.tmp ${gene}.tmp
