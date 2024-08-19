#!/bin/bash
input=$1
module load cd-hit/4.8.1 blast/2.14.0+

#/usr/local/apps/cd-hit/cdhit-4.8.1/psi-cd-hit/psi-cd-hit.pl -i ${input}60 -o ${input}30 -c 0.3
#paste ${input}-clstr-*|cut -f 1-2,4,6 >${input}-clstr.tsv

pre=95
w=5
percent=$(echo ${pre}/100 | bc -l)
cd-hit -i ${input} -o ${input}_${pre} -c ${percent} -n ${w} -d 0 -M 100000 -T 56
python /data/yangy34/projects/useful_scripts/format_cdhit_output.py ${input}_${pre}.clstr ${input}_${pre}.txt
sed -i "s/Cluster /cl${pre}_/" ${input}_${pre}.txt
sort -k4,4 ${input}_${pre}.txt |awk '{print $4"\t"$1}' FS="\t"|sort -k1,1  > ${input}_${pre}_cl.txt

presuffix=_${pre}
for i in 90 85 80 70 60 50 40 30 ;do 
    suffix=${presuffix}-${i}
    percent=$(echo ${i}/100 | bc -l)

    if [ $i -ge 40 ]
    then 
        w=2
    fi

    if [ $i -ge 50 ]
    then 
        w=3
    fi

    if [ $i -ge 60 ]
    then 
        w=4
    fi

    if [ $i -ge 70 ]
    then 
        w=5
    fi

    if [ $i -ge 40 ]
    then 
        cd-hit -i ${input}_${pre} -o ${input}_${i} -c ${percent} -n ${w} -d 0 -M 200000 -T 56
    else
        /usr/local/apps/cd-hit/cdhit-4.8.1/psi-cd-hit/psi-cd-hit.pl -i ${input}_${pre} -o ${input}_${i} -c ${percent}
    fi 

    clstr_rev.pl ${input}${presuffix}.clstr ${input}_${i}.clstr > ${input}${suffix}.clstr
    python /data/yangy34/projects/useful_scripts/format_cdhit_output.py ${input}${suffix}.clstr ${input}${suffix}.txt
    sed -i "s/Cluster /cl${i}_/" ${input}${suffix}.txt
    sort -k4,4 ${input}${suffix}.txt |awk '{print $4"\t"$1}' FS="\t"|sort -k1,1  > ${input}_${i}_cl.txt
    #rm ${input}${pre}.clstr ${input}${pre}  ${input}${presuffix}.clstr ${input}${presuffix}.txt
    pre=${i}
    presuffix=${suffix}
done 


paste *cl.txt|awk '{ print $1; for(i=2;i<=NF;i+=2){print $i};print "\n"}' FS="\t" OFS="" ORS="\t"|sed 's/^\t//;s/\t$//' > ${input}_cluster.txt
rm -rf *.clstr *_cl.txt
