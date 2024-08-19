gid=$1
wd=../data/interproscan
awk -F"\t" '$5=="SSF51126" {print $1"\t"$3"\t"$5}' ${wd}/${gid}/${gid}.*.InterProScan.tsv
