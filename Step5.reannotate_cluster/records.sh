# Step 5. Reannotate biofilm matrix cluster's vps and rbm genes

mkdir -p repr_pff_res_checked
scr=../scripts/curate_Bap1_RbmABC_VpsEF.py
cp ../Step4.annotate_VpsE_VpsF/gid.list .
cat gid.list | parallel -k python ${scr} {} ../Step4.annotate_VpsE_VpsF/repr_pff_res_RbmABC_checked/pffres.{}.annot.gff ../Step4.annotate_VpsE_VpsF/gid_gene_type.RbmABC_Bap1_VpsEF.mapping.pkl ../Step1.annotate_cluster/data/{}.*.gff repr_pff_res_checked/pffres.{}.annot.gff

# Output:
# repr_pff_res_checked/pffres.*.annot.gff
