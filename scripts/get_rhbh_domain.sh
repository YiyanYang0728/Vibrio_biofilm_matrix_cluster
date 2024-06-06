gid=$1
wd=../data/interproscan
awk -F"\t" '$5=="G3DSA:2.160.20.10" {print $1"\t"$3"\t"$5}' ${wd}/${gid}/${gid}.*.InterProScan.tsv