src=/gpfs/gsfs12/users/yangy34/projects/vibrio_YJ_v5/VPS_cluster/scripts/vps_rmb_bap1_count_v3.py
cat /gpfs/gsfs12/users/yangy34/projects/vibrio_paper/Figure1/Fig1A/216_leaves | parallel python ${src} {} repr_pff_res_checked/pffres.{}.annot.gff > 216_gene_ct_v5.tsv
echo -e "gid\tspecies\tvpsA\tvpsB\tvpsD\tvpsE\tvpsF\tvpsI\tvpsJ\tvpsK\tvpsL\tvpsM\tvpsN\tvpsO\trbmA\trbmB\trbmC\trbmD\trbmEF\tBap1\tBap1*" > 216_sp_gene_ct_v5.tsv
paste <(sort -k1,1 /gpfs/gsfs12/users/yangy34/projects/vibrio_paper/Figure1/Fig1A/216_genomes_species.tsv|sed "s/s__//g") <(sort -k1,1 216_gene_ct_v5.tsv)|cut -f1,2,4- >> 216_sp_gene_ct_v5.tsv
/data/yangy34/projects/useful_scripts/binary_table.sh 216_sp_gene_ct_v5.tsv 216_sp_gene_poa_v5.tsv

comp=/gpfs/gsfs12/users/yangy34/projects/vibrio_paper/Figure1/Fig1A/216_sp_gene_poa_v4.tsv
diff 216_sp_gene_poa_v5.tsv ${comp} |grep GCA|awk '{print $2}'|sort -u > diff_gid_poa.list
cat diff_gid_poa.list | parallel /gpfs/gsfs12/users/yangy34/projects/vibrio_paper/Figure1/Fig1A/get_diff_poa_genes.sh {} 216_sp_gene_poa_v5.tsv $comp
# GCA_000695745   rbmB    new:0   old:1
# GCA_003263775   rbmB    new:0   old:1
# GCA_023716585   rbmB    new:0   old:1
# GCA_024347475   rbmB    new:0   old:1
# GCA_921292975   rbmB    new:0   old:1

