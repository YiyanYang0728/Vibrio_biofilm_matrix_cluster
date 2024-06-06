gid=$1
infile=$2

row=`awk '!($0~/#/) && !($0~/Cl_NA/) {print }' ${infile} | wc -l`

if  [ ! ${row} -eq 0 ]
then
    echo ${gid} #awk '!($0~/#/) && !($0~/Cl_NA/) {print }' ${infile}
fi