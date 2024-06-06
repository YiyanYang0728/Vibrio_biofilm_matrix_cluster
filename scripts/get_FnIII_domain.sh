gid=$1
wd=../data/interproscan
awk -F"\t" '$5=="G3DSA:2.60.40.3880"' ${wd}/${gid}/${gid}.*.InterProScan.tsv | cut -f1,3,5,7,8 | sort -u
