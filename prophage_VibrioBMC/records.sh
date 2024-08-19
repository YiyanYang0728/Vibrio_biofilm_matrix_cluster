# classify Vibrio genomes based on completeness of BMC
# 1) have all vps-1, rbmc and vps-2: must have curated vpsEF
# 2) only have part of vps-2: 3 of vpsLMNO
# 3) no BMC clusters: zero vps genes

# first test the representative strains (209+1[GCA_019400765]=210)
file=/gpfs/gsfs12/users/yangy34/projects/vibrio_paper/Figure1/Fig1A/216_sp_gene_poa_v4.tsv
grep -v cholerae ${file} | cat - <(grep GCA_019400765 ${file}) | sed "s/ Clade4//g" > 210_sp_gene_poa.tsv

awk -F"\t" 'NR>1 && $6 > 0 && $7 > 0{print $1}' 210_sp_gene_poa.tsv > complete.gids # must have vpsEF
awk -F"\t" 'NR>1 && $11+$12+$13+$14 >= 3 {print $1}' 210_sp_gene_poa.tsv | grep -vf complete.gids > partial.gids # not in complete list and have 3 of vpsLMNO
grep -vf <(cat complete.gids partial.gids) 210_sp_gene_poa.tsv | grep -wv gid > none.gids # the rest

cat <(awk '{print $1"\tcomplete"}' complete.gids) <(awk '{print $1"\tpartial"}' partial.gids) <(awk '{print $1"\tnone"}' none.gids) > BMC_types.tsv

cut -f1 BMC_types.tsv | parallel -k ../prophage_in_looplessBap1/count_prophages.sh {} > gid_prophage_ct.tsv

echo -e "gid\tbmc_type\tprophage_ct" > BMC_type_prophage_ct.tsv
paste BMC_types.tsv <(cut -f2 gid_prophage_ct.tsv) >> BMC_type_prophage_ct.tsv

echo -e "gid\tbmc_type\tprophage_ct" > BMC_type_vs2_prophage_ct.tsv
awk -F"\t" 'NR==FNR{split($1,x,".");a[x[1]]=$2;next}{print $1"\t"$2"\t"a[$1]}' all_gid_vs2_phagect.tsv BMC_types.tsv >> BMC_type_vs2_prophage_ct.tsv

echo -e "gid\tbmc_type\tprophage_ct" > BMC_type_vibrant_prophage_ct.tsv
awk -F"\t" 'NR==FNR{split($1,x,".");a[x[1]]=$2;next}{print $1"\t"$2"\t"a[$1]}' all_gid_vibrant_phagect.tsv BMC_types.tsv >> BMC_type_vibrant_prophage_ct.tsv
