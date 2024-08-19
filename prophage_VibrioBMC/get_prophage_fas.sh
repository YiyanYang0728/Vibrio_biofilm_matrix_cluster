ctg=$1
start=$2
end=$3
outdir=${4:-$PWD}
gid=`echo ${ctg} | cut -f1-2 -d_`
# echo $gid
cat /gpfs/gsfs12/users/yangy34/projects/vibrio_YJ_v5/VPS_cluster/genomes/fna/${gid}.fna|seqkit fx2tab|fgrep -w ${ctg}|seqkit tab2fx |seqkit subseq -r ${start}:${end} | seqkit fx2tab | awk -v p_ctg=${ctg} -v ns=${start} -v ne=${end} '{print p_ctg"_"ns"-"ne"\t"$2}' | seqkit tab2fx > ${outdir}/${ctg}_${start}_${end}_prophage.fasta
