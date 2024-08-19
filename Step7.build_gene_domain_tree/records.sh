### Step 7. Build gene and domain trees
# 1. Build rbhb (right-handed beta-helix) domain containing gene tree

# Input:
mkdir -p rhbh_tree/
cd rhbh_tree/
# type_mapping.tsv
# type_col.mapping
# genetype2col.mapping
# leaves_phage_asso.tsv

# Process:
# tree
FastTree ../../Step3.annotate_RbmB_RbmA/rhbh_domain_rmdup_trimmed.fasta > rhbh_domain_rmdup.tree

# metadata
gotree labels -i rhbh_domain_rmdup.tree > leaves.list
grep -wf leaves.list ../../Step3.annotate_RbmB_RbmA/rhbh_domain_name_len_label_3.tsv > leaves_info.tsv

awk -F"\t" 'BEGIN{OFS="\t"}NR==FNR{a[$1$2]=$3;next}$3$4 in a{print $0"\t"a[$3$4]}' type_mapping.tsv leaves_info.tsv > leaves_info_type.tsv
join -t$'\t' -1 5 -2 1 <(sort -k5,5 -t$'\t' leaves_info_type.tsv)  <(sort -k1,1 -t$'\t' type_col.mapping) > leaves_info_type_col.tsv

# go to Step8 and get prophage regions first
fgrep -wf leaves.list ../../Step8.find_prophage_regions/leaves_phage_asso.tsv |cut -f1,2 > leaves_phage.tsv
fgrep -wf leaves.list leaves_info_type.tsv |cut -f1,2,5 > leaves_len_genetype.tsv
paste <(sort -k1,1 -t$'\t' leaves_len_genetype.tsv) <(sort -k1,1 -t$'\t' leaves_phage.tsv|cut -f2) > leaves_len_genetype_phage.tsv

# seq len (barplot)
cp bar_tpl.txt bar_annot.txt
awk -F"\t" '{print $1","$2}' leaves_len_genetype_phage.tsv >> bar_annot.txt

# gene type (symbol at leaves)
cp symbol_tpl.txt symbol_annot.txt
awk -F"\t" 'NR==FNR{a[$1]=$2;next}{print $1",2,5,"a[$3]",1,1,"$3}' genetype2col.mapping leaves_len_genetype_phage.tsv >> symbol_annot.txt

# phage closeness (color strip)
cp strip_tpl.txt strip_annot.txt
awk -F"\t" '$4=="In_prophage"{print $1",#2C3279"}$4=="Not_in_prophage"{print $1",#B8C4D5"}$4=="Same_contig"{print $1",#03B3A9"}' leaves_len_genetype_phage.tsv >> strip_annot.txt

zip -r rhbh_domain_rmdup.tree.zip rhbh_domain_rmdup.tree bar_annot.txt symbol_annot.txt strip_annot.txt
upload_tree rhbh_domain_rmdup.tree.zip
# https://itol.embl.de/tree/128231274380931723775282
# copy 143_phage.genes names from iTOL tree

# blast genes against INPHARED
cat 143_phage.genes | parallel /data/yangy34/projects/vibrio_YJ_v5/VPS_cluster/genomes/get_faa_seq.sh {} > phage.genes.faa
blastp -num_threads 32 -query phage.genes.faa -subject /gpfs/gsfs12/users/yangy34/database/inphared/2Aug2024_vConTACT2_proteins.faa -outfmt 6 -out inphared.blast.out

cd ..

# Output:
# rhbh_domain_rmdup.tree.zip

# 2. Build beta-propeller domain tree
# Beta-propeller domain tree

# Input:
mkdir -p b-propeller_tree
cd b-propeller_tree/
cat ../../Step2.annotate_RbmC_Bap1/all.ids ../../Step2.annotate_RbmC_Bap1/RbmC_w_b-helix.ids > gene_ids.list
# get beta-prop domain fasta
seqkit fx2tab ../../Step2.annotate_RbmC_Bap1/all_genes.rmdup.aln.faa | awk '{print $1"\t"substr($NF,244,552-244+1)substr($NF,766,968-766+1)}' | seqkit tab2fx -w 0 > Beta_propeller.1.fasta
# get RbmC with beta-helix domain proteins' beta-prop domain fasta
seqkit fx2tab ../../Step2.annotate_RbmC_Bap1/all.faa | fgrep -wf ../../Step2.annotate_RbmC_Bap1/RbmC_w_b-helix.ids |awk '{print $1"\t"substr($NF,286,1000)}' | seqkit tab2fx -w 0 > RbmC_w_b-helix_beta_propeller.fasta

# Process:
# align 997+5 domains
cat <(sed "s/-//g" Beta_propeller.1.fasta) RbmC_w_b-helix_beta_propeller.fasta > Beta_propeller.unaligned.fasta
mafft --thread 32 --maxiterate 1000 --localpair Beta_propeller.unaligned.fasta > Beta_propeller.3.fasta

seqkit fx2tab Beta_propeller.3.fasta | awk -F"\t" '{print $1"\t"substr($2,14,length($2)-(14+7)+1)}' | seqkit tab2fx -w 0 > Beta_propeller.fasta

# tree
FastTree Beta_propeller.fasta > Beta_propeller.tree

# species info
grep -v GCA ../../Step2.annotate_RbmC_Bap1/itol_tree/taxon_annot.txt > taxon_tpl.txt
cp taxon_tpl.txt taxon_annot.txt
grep -wf gene_ids.list ../../Step2.annotate_RbmC_Bap1/itol_tree/taxon_annot.txt >> taxon_annot.txt
echo -e "GCA_003544875.1_03596,#A39BA8,s__Vibrio alfacsensis\nGCA_019670025.1_03371,#A39BA8,s__Vibrio alfacsensis\nGCA_019670485.1_03306,#A39BA8,s__Vibrio alfacsensis\nGCA_905175395.2_03458,#A39BA8,s__Vibrio alfacsensis\nGCA_002608565.1_00214,#EDF5FC,s__Vibrio sp002608565" >> taxon_annot.txt

# structure info
grep -v GCA ../../Step2.annotate_RbmC_Bap1/itol_tree/struct_type_annot.2.txt > struct_tpl.txt
cp struct_tpl.txt struct_type_annot.txt
grep -wf gene_ids.list ../../Step2.annotate_RbmC_Bap1/itol_tree/struct_type_annot.2.txt >> struct_type_annot.txt
echo  -e "GCA_003544875.1_03596,#69146B,RbmC_with_b-helix\nGCA_019670025.1_03371,#69146B,RbmC_with_b-helix\nGCA_019670485.1_03306,#69146B,RbmC_with_b-helix\nGCA_905175395.2_03458,#69146B,RbmC_with_b-helix\nGCA_002608565.1_00214,#69146B,RbmC_with_b-helix" >> struct_type_annot.txt

zip -r Beta_propeller.tree.zip Beta_propeller.tree struct_type_annot.txt taxon_annot.txt
upload_tree Beta_propeller.tree.zip

cd ..

# Output:
# Beta_propeller.tree.zip

# 3. Build beta-prism tree

# Input:
mkdir -p b-prism_tree
cd b-prism_tree/

cp ../../Step2.annotate_RbmC_Bap1/RbmC.ids RbmC.list
cp ../../Step2.annotate_RbmC_Bap1/Bap1.ids Bap1.list

seqkit fx2tab ../../Step2.annotate_RbmC_Bap1/all_genes.rmdup.aln.faa |grep -wf <(cat RbmC.list Bap1.list) | awk '{print $1"_prismC1\t"substr($NF,553,659-553+1)substr($NF,717,765-717+1)}' | sed "s/-//g" | awk 'length($2)>110' | seqkit tab2fx -w 0 > Beta_prisimC1.fasta
seqkit fx2tab ../../Step2.annotate_RbmC_Bap1/all_genes.rmdup.aln.faa |grep -wf RbmC.list | awk '{print $1"_prismC2\t"substr($NF,969,1160-969+1)}'  | sed "s/-//g" | awk 'length($2)>110' | seqkit tab2fx -w 0 > Beta_prisimC2.fasta

# Process:
cat Beta_prisimC1.fasta Beta_prisimC2.fasta > Beta_prisims.unaligned.fasta
mafft --thread 32 --maxiterate 1000 --localpair Beta_prisims.unaligned.fasta > Beta_prisims.fasta

#  tree
FastTree Beta_prisims.fasta > Beta_prisims.tree

# species info
seqkit fx2tab Beta_prisims.fasta|cut -f1 > gene_ids.list
sed "s/_prism/\tprism/g" gene_ids.list > leaves_prism_type.tsv
grep -v GCA ../../Step2.annotate_RbmC_Bap1/itol_tree/taxon_annot.txt > taxon_tpl.txt
cp taxon_tpl.txt taxon_annot.txt
awk -F"\t" 'NR==FNR{if($0~/^GCA/){split($0,a,",");b[a[1]]=a[2]","a[3]};next}$1 in b{print $1"_"$2","b[$1]}' ../b-propeller_tree/taxon_annot.txt leaves_prism_type.tsv >> taxon_annot.txt

# structure info
grep -v GCA ../../Step2.annotate_RbmC_Bap1/itol_tree/struct_type_annot.2.txt > struct_tpl.txt
cp struct_tpl.txt struct_type_annot.txt
awk -F"\t" 'NR==FNR{if($0~/^GCA/){split($0,a,",");b[a[1]]=a[2]","a[3]};next}$1 in b{print $1"_"$2","b[$1]}' ../b-propeller_tree/struct_type_annot.txt leaves_prism_type.tsv >> struct_type_annot.txt

zip -r Beta_prisims.tree.zip Beta_prisims.tree struct_type_annot.txt taxon_annot.txt
upload_tree Beta_prisims.tree.zip

cd ..

# Output:
# Beta_prisims.tree.zip