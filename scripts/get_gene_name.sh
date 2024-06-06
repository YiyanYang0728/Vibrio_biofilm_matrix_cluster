gene=$1
gid=$(echo ${gene}|cut -f1 -d.)
gene_name=`fgrep -w ${gene} -m 1 ../Step1.annotate_cluster/data/${gid}.*.gff|awk -F"\t" '{split($NF, a, ";");for(i=1;i<=length(a);i=i+1){if(a[i]~/product/){split(a[i],b,"=");print b[2]}}}'`
echo -e "${gene}\t${gene_name}"
