gene=$1
gid=`echo ${gene}|cut -f1 -d.`
wd=../data/interproscan
awk -F"\t" '($5=="PF13440" || $5=="PF14667")' <(fgrep ${gene} ${wd}/${gid}/${gid}.*.InterProScan.tsv) | cut -f1,3,5,7,8 | sort -u
