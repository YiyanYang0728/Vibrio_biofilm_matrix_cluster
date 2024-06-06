### Step 4. Annotate VpsE and VpsF proteins
mkdir -p repr_pff_res_RbmABC_checked
scr=../scripts/curate_Bap1_RbmABC_VpsEF.py
cp ../Step3.annotate_RbmB_RbmA/gid.list .
cat gid.list | parallel -k python ${scr} {} ../Step3.annotate_RbmB_RbmA/repr_pff_res_RbmBC_checked/pffres.{}.annot.gff ../Step3.annotate_RbmB_RbmA/gid_gene_type.RbmABC_Bap1.mapping.pkl ../Step1.annotate_cluster/data/{}.*.gff repr_pff_res_RbmABC_checked/pffres.{}.annot.gff

cat gid.list | parallel -k "grep 'Name=vpsE' repr_pff_res_RbmABC_checked/pffres.{}.annot.gff|grep -v "ClusterID=Cl_NA" |grep -oP 'ID=\KGCA_[0-9]*\.[0-9]*_[0-9]*'" > vpsE_ids.list
cat gid.list | parallel -k "grep 'Name=vpsF' repr_pff_res_RbmABC_checked/pffres.{}.annot.gff|grep -v "ClusterID=Cl_NA" |grep -oP 'ID=\KGCA_[0-9]*\.[0-9]*_[0-9]*'" > vpsF_ids.list

cat vpsE_ids.list | parallel -k ../scripts/get_PF13440_PF14667_domain.sh {} > vpsE_ids_domain.list
cat vpsF_ids.list | parallel -k ../scripts/get_NF038256_domain.sh {} > vpsF_ids_domain.list

# # 1.1 get vpsE standard list
cut -f1 vpsE_ids_domain.list|sort|uniq -c|awk '$1==2{split($NF, a, ".");print a[1]"\t"$NF}' > vpsE_std.list # 1749, $2==2 means to have both
# # 1.2 get hmm failed-to-recognize vpsE ==> write to vpsE_suppl_1.list
# GCA_000338875.1_03714
# GCA_003311865.1_03558
# GCA_006802925.1_00025
# GCA_013414565.1_00448
# GCA_006802905.1_00928
# # 1.3 multi-gene vpsE ==> write to vpsE_suppl_2.list
# GCA_000168935.3_00011 + GCA_000168935.3_00012 (13+14)
# GCA_003056705.1_02300 + GCA_003056705.1_02301 (13+14)
# GCA_000961975.1_02072 + GCA_000961975.1_02073 + GCA_000961975.1_02074 (13+13+14)
# GCA_001637545.1_02254 + GCA_001637545.1_02256 (13+14)
# GCA_001637555.1_03503 + GCA_001637555.1_02683 (13+14)

awk '{split($NF, a, ".");print a[1]"\t"$NF}' vpsE_suppl_1.list > vpsE_suppl_1.1.list
awk '{split($NF, a, ".");print a[1]"\t"$NF}' vpsE_suppl_2.list > vpsE_suppl_2.1.list
cat vpsE_std.list vpsE_suppl_1.1.list vpsE_suppl_2.1.list > vpsE_all.list

# # 2.1 get vpsF standard list
cut -f1 vpsF_ids_domain.list |sort -u| awk '{split($NF, a, ".");print a[1]"\t"$NF}' > vpsF_std.list # 1765
# # 2.2 get hmm failed-to-recognize vpsF/multi-gene vpsF ==> write to vpsF_suppl_1.list
# GCA_000237745.2_01774 + GCA_000237745.2_01773
# GCA_002196325.1_05261
# GCA_006802795.1_02354
# # 2.3 get domain failed vpsF ==> write to vpsF_suppl_2.list
# GCA_011212705.1_03731

awk '{split($NF, a, ".");print a[1]"\t"$NF}' vpsF_suppl_1.list > vpsF_suppl_1.1.list
awk '{split($NF, a, ".");print a[1]"\t"$NF}' vpsF_suppl_2.list > vpsF_suppl_2.1.list
cat vpsF_std.list vpsF_suppl_1.1.list vpsF_suppl_2.1.list > vpsF_all.list

# create mapping
cut -f2 vpsE_all.list | tr "," "\n"|awk -F"\t" '{split($1,a,"_"); print a[1]"_"a[2]"\t"$1"\tVpsE"}' > vpsE.mapping
cut -f2 vpsF_all.list | tr "," "\n"|awk -F"\t" '{split($1,a,"_"); print a[1]"_"a[2]"\t"$1"\tVpsF"}' > vpsF.mapping
cat vpsE.mapping vpsF.mapping ../Step3.annotate_RbmB_RbmA/gid_gene_type.RbmABC_Bap1.mapping > gid_gene_type.RbmABC_Bap1_VpsEF.mapping
python ../scripts/pickle_file.py gid_gene_type.RbmABC_Bap1_VpsEF.mapping

# Output:
# gid_gene_type.RbmABC_Bap1_VpsEF.mapping
