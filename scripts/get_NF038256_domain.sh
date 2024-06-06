gene=$1
gid=`echo ${gene}|cut -f1 -d.`
wd=../data/interproscan
awk -F"\t" '$5=="NF038256"' <(fgrep ${gene} ${wd}/${gid}/${gid}.*.InterProScan.tsv) | cut -f1,3,5,7,8 | sort -u