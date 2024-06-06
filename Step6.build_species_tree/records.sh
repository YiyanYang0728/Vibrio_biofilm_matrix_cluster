### Step 6. Pick repr genomes for Vibrio species and build species tree
# 1. select 210 species representative genomes

# Input:
# gid.list
cp ../Step3.annotate_RbmB_RbmA/gid.list .
# gid_sp.tsv
grep -wf gid.list /gpfs/gsfs12/users/yangy34/projects/vibrio_YJ_v5/GTDB_r214_taxon.tsv|awk -F"\t" '{split($3,a,";");split($2,b,".");print b[1]"\t"a[7]}'| sed "s/s__//g" > gid_sp.tsv
# Manually create files for the selected genomes as below:
# 22_repr_genomes.tsv
# 7_Vc_repr_genomes.tsv

# Process:
# 187 species (4111 genomes) have no repr genomes
cat <(cut -f2 22_repr_genomes.tsv) <(echo "Vibrio cholerae")|grep -vwf - gid_sp.tsv > gid_sp_rest.tsv

# for genomes in gid_sp_rest.tsv (187 species), only keep genomes have clustered VPS genes
src=../scripts/genome_w_cluster.sh
cut -f1 gid_sp_rest.tsv | parallel ${src} {} ../Step5.reannotate_cluster/repr_pff_res_checked/pffres.{}.annot.gff > gid_sp_rest.VPS.gids # 3315
grep -wf gid_sp_rest.VPS.gids gid_sp_rest.tsv > gid_sp_rest_VPS.tsv # 3315 genomes in 114 species
cut -f2 gid_sp_rest_VPS.tsv|sort -u > gid_sp_rest_VPS.sp # 114 species that have at least one genome with VPS cluster
grep -vwf gid_sp_rest_VPS.sp gid_sp_rest.tsv > gid_sp_rest_noVPS.tsv # 157 genomes in 73 species have no VPS
grep -wf <(cut -f2 gid_sp_rest_noVPS.tsv) repr_gid2taxon.tsv > gid_sp_rest_noVPS_reprgid_sp.tsv # 73 species with GTDB repr genomes
grep -wf gid_sp_rest_VPS.sp gid_sp_rest.tsv | cut -f2 |sort |uniq -d|grep -wf - gid_sp_rest.tsv > gid_sp_rest_VPS_multigid.tsv # 3916 genomes of 76 species
grep -wf gid_sp_rest_VPS.sp gid_sp_rest.tsv | cut -f2 |sort |uniq -u|grep -wf - gid_sp_rest.tsv > gid_sp_rest_VPS_singlegid.tsv # 38 genomes of 38 species
# 114 species have potential VPS clusters: 76 sp with >=2 genomes, 38 sp with 1 genome

# Pick repr genome for these 61 species (with >=2 genomes)
# for each species, pick repr genomes 1) have the least contigs (most complete); 2) have the biggest VPS cluster; 3) better to have rbm/bap1 genes
# ==> create a table: contigs #; vps-I cluter Very Important Genes # (ABDEFIJK); vps-II cluter Very Important Genes # (LMNO); rbm genes #
cut -f1 gid_sp_rest_VPS_multigid.tsv | parallel -k ../scripts/count_contig.sh {} > gid_sp_rest_VPS_multigid.contig_ct

src=../scripts/vps_rmb_count.py
cut -f1 gid_sp_rest_VPS_multigid.tsv | parallel python ${src} {} ../Step5.reannotate_cluster/repr_pff_res_checked/pffres.{}.annot.gff > gid_sp_rest_VPS_multigid.gene_ct

echo -e "gid\tspecies\tcontig_ct\tvps-I_VIG_ct\tvps-II_VIG_ct\trbm_genes" > gid_sp_rest_VPS_multigid_info.tsv
paste <(sort -k1,1 gid_sp_rest_VPS_multigid.tsv) <(sort -k1,1 gid_sp_rest_VPS_multigid.contig_ct |cut -f2) <(sort -k1,1 gid_sp_rest_VPS_multigid.gene_ct|cut -f2-) >> gid_sp_rest_VPS_multigid_info.tsv

tail -n+2 gid_sp_rest_VPS_multigid_info.tsv|sort -t$'\t' -k2,2 -k3,3n -k5,5nr -k4,4nr|awk -F"\t" '!seen[$2]++' > gid_sp_rest_VPS_repr_gid_info.tsv

# 73 species repr genomes with no VPS: cut -f1,2 gid_sp_rest_noVPS_reprgid_sp.tsv
# 76 species repr genomes: cut -f1,2 gid_sp_rest_VPS_repr_gid_info.tsv
# 38 species repr genomes: cut -f1,2 gid_sp_rest_VPS_singlegid.tsv
# 22 species repr genomes: cut -f1,2 22_repr_genomes.tsv
# 73+76+38+22+1=210
cat <(cut -f1,2 gid_sp_rest_noVPS_reprgid_sp.tsv) <(cut -f1,2 gid_sp_rest_VPS_repr_gid_info.tsv) \
 <(cut -f1,2 gid_sp_rest_VPS_singlegid.tsv) <(cut -f1,2 22_repr_genomes.tsv) <(echo -e "GCA_019400765\tVibrio cholerae") > 210_repr_gid_sp.tsv
 
# 7 Vc subspecies repr genomes: cut -f1,2 7_Vc_repr_genomes.tsv
# A total of 216 (sub)species repr genomes
# rest species have no clusters

# 2. build species trees
cat <(grep -v cholerae 210_repr_gid_sp.tsv) 7_Vc_repr_genomes.tsv > 216_repr_gid_sp.tsv
mkdir -p gff
cut -f1 216_repr_gid_sp.tsv | parallel -k cp ../Step1.annotate_cluster/data/{}.*.gff gff/
PIRATE -i gff/  -o pirate_out/ -a -t 32 -k "--diamond"
FastTree -gtr -nt pirate_out/core_alignment.fasta > pirate_out/core_alignment.tree

# get 29 species + outgroup substree
cat 22_repr_genomes.tsv 7_Vc_repr_genomes.tsv <(grep Vibrio_A 216_repr_gid_sp.tsv) |cut -f1 > 30_leaves
gotree labels -i pirate_out/core_alignment.tree > pirate_out/core_alignment.leaves
gotree prune -r -f <(grep -f 30_leaves pirate_out/core_alignment.leaves) -i pirate_out/core_alignment.tree -o 30_core_alignment.sub.tree

# Output:
# 216_repr_gid_sp.tsv
# pirate_out/core_alignment.tree
# 30_core_alignment.sub.tree