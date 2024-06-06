gid=$1
contig_ct=`grep "##sequence-region" ../Step1.annotate_cluster/data/${gid}.*.gff|wc -l`
echo -e "${gid}\t${contig_ct}"
