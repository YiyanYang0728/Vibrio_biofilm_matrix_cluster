geneId=$1
gid=`echo ${geneId}|cut -f1 -d. `
grep -w ${geneId} ../Step1.annotate_cluster/data/${gid}.*.gff |awk -v geneid=${geneId} '{print geneid"\t"($5-$4+1)/3-1}'