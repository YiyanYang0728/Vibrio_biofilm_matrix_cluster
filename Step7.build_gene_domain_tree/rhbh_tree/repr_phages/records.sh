# get gene and phage regions
grep -wf ../143_phage.genes ../../../Step8.find_prophage_regions/leaves_phage_asso.tsv|\
awk '{print $1"\t"$5"\t"$NF}' | grep -v contig > gene_loc_phage.tsv

# get phage region fasta
get_phage=/data/yangy34/projects/vibrio_paper/prophage_VibrioBMC/get_prophage_fas.sh
cut -f1,3 gene_loc_phage.tsv|awk '{split($2,a,":");split(a[2],b,"-");print $1"\t"a[1]"\t"b[1]"\t"b[2]}' > gene_phage_region.tsv
mkdir -p phage_fas 
cut -f2- gene_phage_region.tsv |parallel --colsep '\t' ${get_phage} {1} {2} {3} phage_fas

# cluster phage genomes and get repr
cat phage_fas/*.fasta > phages.fas
/data/yangy34/softwares/MGV_cluster_phage/cluster_phage_ident70.sh phages.fas

# 13 repr phages
# wc -l clusters_70.tsv
cut -f1 clusters_70.tsv > repr_phages.list
# add 2 V. metoecus and V. sp000176715 genomes
echo "GCA_001402675.1_4_217771-307261" >> repr_phages.list
echo "GCA_000176715.1_9_60706-133991" >> repr_phages.list

mkdir -p phage_pharokka
cat repr_phages.list | parallel "cp -r /gpfs/gsfs12/users/yangy34/projects/Vibrio_phages/run_vs2/phage_pharokka/{} phage_pharokka"

mkdir -p clinker/gff
for id in `cat repr_phages.list`
do
    cp phage_pharokka/${id}/${id}.gff clinker/gff/
    id2=`echo $id|sed "s/-/_/g"`
    cp phage_fas/${id2}_prophage.fasta clinker/gff/${id}.fna
done

# create gene2grp.csv
cp /gpfs/gsfs12/users/yangy34/projects/vibrio_YJ_v5/VPS_cluster/genomes/RbmB_genes/Vibrio_rhbh_evol/prophage_annot/find_target_in_prophage_v2.sh find_target_in_prophage.sh
cut -f1 -d. repr_phages.list|grep -wf - gene_phage_region.tsv|cut -f1| parallel -k ./find_target_in_prophage.sh {} | sort -u -k1,1 > gene2alias.mapping
cat repr_phages.list| parallel "tail -n+2 phage_pharokka/{}/{}_out.tsv | cut -f1,16" | cut -f1 -d"," | tr "\t" "," > gene2ctg.tsv
awk -F"," 'NR==FNR{a[$1];next}$1 in a{print $1",Rhbh"}!($1 in a){print $1","$2}' <(cut -f2 gene2alias.mapping) gene2ctg.tsv > gene2grp.csv
# create grp2col.csv
cp /data/yangy34/projects/vibrio_YJ_v5/VPS_cluster/genomes/RbmB_genes/Vibrio_rhbh_evol/prophage_annot/grp2col.csv .

cd clinker
clinker gff/*.gff --no_align -f -gf ../gene2grp.csv -cm ../grp2col.csv -p 15_prophage_plot.html