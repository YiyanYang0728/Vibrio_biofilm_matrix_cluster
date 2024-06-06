gene=$1
gid=$(echo ${gene}|cut -f1 -d.)
wd=../Step1.annotate_cluster/data
seqkit fx2tab ${wd}/${gid}*.faa | fgrep ${gene}|awk '{print $1"\t"$NF}'| seqkit tab2fx -w 0
