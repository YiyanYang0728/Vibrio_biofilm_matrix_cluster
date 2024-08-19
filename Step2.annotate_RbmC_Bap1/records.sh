### Step 2. Annotate RbmC adn Bap1 proteins
# Input:
# Bap1_RbmC.fa
# all.faa
# all.ffn
# all.faa.* [blast db file]
# gid2gtdb.tsv 
# gtdb_taxon.mapping

# Process:
### 1) get putative RbmC and Bap1 genes
blastp -query Bap1_RbmC.fa -db all.faa -outfmt 6 -num_threads 32 -max_target_seqs 26240000 > query_all_Bap1_RbmC.out
awk -F"\t" '$3>40 && $4>200 && $NF>250' query_all_Bap1_RbmC.out > query_all_Bap1_RbmC_filtered.out
cut -f2 query_all_Bap1_RbmC_filtered.out | sort -u > query_Bap1_RbmC_gene.list
# remove GCA_002933855.1_04065 and GCA_008084795.1_02439 and RbmC_w_b-helix in the list
grep -v GCA_002933855.1_04065 query_Bap1_RbmC_gene.list | grep -v GCA_008084795.1_02439 | grep -vf RbmC_w_b-helix.ids > query_Bap1_RbmC_gene.filtered.list
seqkit grep -f query_Bap1_RbmC_gene.filtered.list all.faa -j 32 -o all_genes.faa # 4059
seqkit grep -f query_Bap1_RbmC_gene.filtered.list all.ffn -j 32 -o all_genes.ffn # 4059
seqkit rmdup -s all_genes.ffn -d all_genes.dup.ffn -D all_genes.dup.ffn.ids > all_genes.rmdup.ffn # 997 seqs
seqkit grep -n -f <(grep ">" all_genes.rmdup.ffn|sed "s/>//g") all_genes.faa -j 32 -o all_genes.rmdup.faa # 997 seqs

### 2) create MSA
awk '/^>/{print s? s"\n"$0:$0;s="";next}{s=s sprintf("%s",$0)}END{if(s)print s}' all_genes.rmdup.ffn|awk 'NR%2==1{name=$1}NR%2==0{print name"\t"length($0)}'|sed "s/>//g" > all_genes_id_len.rmdup.tsv
# Bap1
cut -f1 all_genes_id_len.rmdup.tsv |grep -wf - query_all_Bap1_RbmC_filtered.out |awk '$1=="Bap1" && $4>650 && $4<700 && $3>80 && $NF>900{print $2}'|sort|uniq > high_quality_id.list
# RbmC
cut -f1 all_genes_id_len.rmdup.tsv |grep -wf - query_all_Bap1_RbmC_filtered.out |awk '$1=="RbmC" && $4>950 && $4<1000 && $3>80 && $NF>900{print $2}'|sort|uniq >> high_quality_id.list
cut -f1 all_genes_id_len.rmdup.tsv |sort high_quality_id.list -|uniq -u > low_quality_id.list

seqkit grep -f high_quality_id.list all_genes.rmdup.faa -j 32 -o high_quality_id.faa # 672
seqkit grep -f low_quality_id.list all_genes.rmdup.faa -j 32 -o low_quality_id.faa # 325

mafft --thread 32 --maxiterate 1000 --localpair high_quality_id.faa > high_quality_id_aln.faa

# 1-2 use mafft --add to add low-quality sequences
mafft --add low_quality_id.faa --thread 32 --maxiterate 1000 --localpair high_quality_id_aln.faa > all_genes.rmdup.aln.faa

# test alignment by checking Bap1: GCA_003347715.1_00465 and RbmC: GCA_004801065.1_02413
awk '/^>/{print s? s"\n"$0:$0;s="";next}{s=s sprintf("%s",$0)}END{if(s)print s}' all_genes.rmdup.aln.faa|grep -P "GCA_003347715.1_00465|GCA_004801065.1_02413" -A 1 --no-group-separator > test_aln.fa
# test if 57aa loop is aligned well
grep YLGLEWKTKTVPYLGVEWRTKTVSYWFFGWHTKQVAYLAPVWKEKTIPYAVPVTLSK test_aln.fa -C 2

python ../scripts/fas2clustal.py all_genes.rmdup.aln.faa all_genes.rmdup.aln
/data/yangy34/softwares/pal2nal.v14/pal2nal.pl all_genes.rmdup.aln all_genes.rmdup.ffn -output fasta -codontable 11 > all_genes.rmdup.codon.aln.fa

### 3) build tree
# create dna12_3.partition.txt like this:
# DNA, p1=1-3480\3,2-3480\3
# DNA, p2=3-3480\3

ml raxml
mkdir -p run_raxml
cd run_raxml
cp ../dna12_3.partition.txt .
cp ../all_genes.rmdup.codon.aln.fa .
raxmlHPC-PTHREADS-SSE3 -T 56 -m GTRGAMMA -p 12345 -q dna12_3.partition.txt -s all_genes.rmdup.codon.aln.fa -n raxml.v2.out
# final result tree:
# run_raxml/RAxML_bestTree.raxml.v2.out
cd ..

# get 57aa loop status
echo ">test" > test.fa
grep  TATCTAGGATTAGAGTGGAAAACTAAAACGGTTCCTTATCTAGGTGTTGAGTGGCGTACCAAAACCGTCTCTTACTGGTTCTTTGGCTGGCACACTAAACAAGTGGCTTATTTAGCGCCAGTATGGAAAGAGAAAACCATCCCTTATGCCGTTCCTGTGACACTGTCGAAA <(seqkit fx2tab all_genes.rmdup.codon.aln.fa|seqkit tab2fx -w 0) -m 1 >> test.fa
seqkit locate -p TATCTAGGATTAGAGTGGAAAACTAAAACGGTTCCTTATCTAGGTGTTGAGTGGCGTACCAAAACCGTCTCTTACTGGTTCTTTGGCTGGCACACTAAACAAGTGGCTTATTTAGCGCCAGTATGGAAAGAGAAAACCATCCCTTATGCCGTTCCTGTGACACTGTCGAAA test.fa|cut -f5,6
# 57aa/171nt is located at
# start   end
# 1978    2148

goalign subsites -i test.fa {1977..2147} | seqkit translate -T 11 # checked: the same as "YLGLEWKTKTVPYLGVEWRTKTVSYWFFGWHTKQVAYLAPVWKEKTIPYAVPVTLSK"
goalign subsites -i <(seqkit fx2tab all_genes.rmdup.codon.aln.fa|seqkit tab2fx -w 0) {1977..2147} > 171nt.fa
seqkit translate -T 11 -j 32 -o 171nt_translated.fa 171nt.fa
awk '/^>/{print s? s"\n"$0:$0;s="";next}{s=s sprintf("%s",$0)}END{if(s)print s}' 171nt.fa |awk 'NR%2==1{name=$1}NR%2==0{print name"\t"$0}'|sed "s/>//g" > id_171nt.tsv
awk '/^>/{print s? s"\n"$0:$0;s="";next}{s=s sprintf("%s",$0)}END{if(s)print s}' 171nt_translated.fa |awk 'NR%2==1{name=$1}NR%2==0{print name"\t"$0}'|sed "s/>//g" > id_171nt_trans.tsv
cut -f2 id_171nt.tsv|sed "s/-//g"|awk '{print length($0)}' | paste id_171nt.tsv <(cut -f2 id_171nt_trans.tsv) - > id_171nt_validlen.tsv
awk -F"\t" 'BEGIN{OFS="\t"}$4==0{label="empty"}$4>0 && $4<171{label="not_full"}$4==171{label="full"}{print $0"\t"label}' id_171nt_validlen.tsv > id_171nt_validlen_label.tsv

awk 'BEGIN{OFS="\t"}{a=substr($1,1,13);print a}' id_171nt_validlen_label.tsv | grep -wf - <(cut -f2-3 gid2gtdb.tsv) | sed "s/\.[0-9]*//g" | awk '{print $2"\t"$1}' > gtdb_gid.mapping
join -t$'\t' -1 2 -2 1 <(sort -k2,2 gtdb_gid.mapping) <(sort -k1,1 gtdb_taxon.mapping) > gtdb_gid_taxon_col.tsv
awk 'BEGIN{OFS="\t"}{a=substr($1,1,13);print a,$0}' id_171nt_validlen_label.tsv > gid_id_171nt_validlen_label.tsv
join -t$'\t' -1 1 -2 2 <(sort -k1,1 gid_id_171nt_validlen_label.tsv) <(sort -k2,2 gtdb_gid_taxon_col.tsv) > gid_id_171nt_validlen_label_gtdb_species_col.tsv
awk -F"\t" 'BEGIN{OFS="\t"}{print $2,$(NF-1),$NF}' gid_id_171nt_validlen_label_gtdb_species_col.tsv > id_taxon_col.tsv

mkdir -p itol_tree
cd itol_tree/
cp ../dataset_symbols_template.txt node_annot.txt
awk -F"\t" 'BEGIN{OFS=","}$5=="full"{print $1,"2",10,"#F46036",1,0}$5=="not_full"{print $1,"2",10,"#1B998B",1,0}$5=="empty"{print $1,"2",10,"#C5D86D",1,0}' ../id_171nt_validlen_label.tsv >> node_annot.txt

cp ../dataset_color_strip_template.txt taxon_annot.txt
awk -F"\t" 'BEGIN{OFS=","}{print $1,$3,$2}' ../id_taxon_col.tsv >> taxon_annot.txt

cp ../run_raxml/RAxML_bestTree.raxml.out gene.tree

zip -r gene.tree.zip gene.tree node_annot.txt taxon_annot.txt
upload_tree gene.tree.zip # 997 leaves
# itol link: https://itol.embl.de/tree/128231272458421717441471

### 4) get Bap1, loop-less Bap1 and RbmC labels
# from tree manually get all.ids and Bap1.ids
grep -vf Bap1.ids all.ids > RbmC.ids
cut -f1,4,5 id_171nt_validlen_label.tsv > genes_57aa_poa.tsv
grep -wf Bap1.ids genes_57aa_poa.tsv | awk '($2/171 < 0.12){print $1}' > Bap1_loopless.ids
# check that 21 (GCA_000176455.1_03193),45 (GCA_013155105.1_06224),69 (GCA_001860485.1_01253) are all **partial genes** with only a b-propellor domain as classified as Bap1
grep -vf Bap1_loopless.ids Bap1.ids > rest_Bap1.ids
### summary:
#    107 Bap1_loopless.ids
#    514 RbmC.ids
#      5 RbmC_w_b-helix.ids
#    376 rest_Bap1.ids

# merge gene types and set GCA_002933855.1_04065 as RbmC_fake
cat <(awk '{print $0"\tRbmC"}' RbmC.ids) <(awk '{print $0"\tBap1"}' rest_Bap1.ids) <(awk '{print $0"\tBap1_wo_57aa"}' Bap1_loopless.ids) <(awk '{print $0"\tRbmC_w_b-helix"}' RbmC_w_b-helix.ids) <(echo -e "GCA_002933855.1_04065\tRbmC_fake") > gene_types.1.tsv
# first rough typing
cut -f2 gene_types.1.tsv|sort |uniq -c
# 376 Bap1
# 107 Bap1_wo_57aa
# 514 RbmC
#   1 RbmC_fake
#   5 RbmC_w_b-helix

### 5) check by Geneious to get domain boundaries
# all_genes.rmdup.aln_Annotations.tsv

### 6) type genes
python ../scripts/gene_typing.py all_genes.rmdup.aln.faa <(cut -f1 gene_types.1.tsv) gene_types.1.tsv gene_types.2.tsv
python ../scripts/assign_domain_presabs.py all_genes.rmdup.aln.faa > all_genes.domain_presabs.tsv # for FigS3
# second refined typing
cut -f2 gene_types.2.tsv|sort |uniq -c
#     358 Bap1_standard
#      18 Bap1_truncated
#      11 Bap1_wo_57aa_truncated
#      96 Bap1_w/o_loop
#       1 RbmC_fake
#     428 RbmC_standard
#      44 RbmC_truncated
#       5 RbmC_with_b-helix
#      20 RbmC_with_partial_M1M2
#      22 RbmC_w/o_M1M2

cut -f2- all_genes.dup.ffn.ids|awk '{split($0,a,", ");for(i=2;i<=length(a);i=i+1){print a[1]"\t"a[i]}}' > repr2dup.mapping
awk -F"\t" 'NR==FNR{a[$1]=$2}NR!=FNR{if($1 in a){print $2"\t"a[$1]}}' gene_types.2.tsv repr2dup.mapping > gene_types.rest.tsv
cat gene_types.2.tsv gene_types.rest.tsv > gene_types.all.tsv # 4065
cut -f2 gene_types.all.tsv|sort |uniq -c
#    1784 Bap1_standard
#      29 Bap1_truncated
#      11 Bap1_wo_57aa_truncated
#     238 Bap1_w/o_loop
#       1 RbmC_fake
#    1855 RbmC_standard
#      85 RbmC_truncated
#       5 RbmC_with_b-helix
#      26 RbmC_with_partial_M1M2
#      31 RbmC_w/o_M1M2

paste <(cut -f1 gene_types.all.tsv|cut -f1-2 -d_) gene_types.all.tsv > gid_gene_type.RbmC_Bap1.mapping
python ../scripts/pickle_file.py gid_gene_type.RbmC_Bap1.mapping

cut -f3 gid_gene_type.RbmC_Bap1.mapping|sort|uniq -c
#    1784 Bap1_standard
#      29 Bap1_truncated
#      11 Bap1_wo_57aa_truncated
#     238 Bap1_w/o_loop
#       1 RbmC_fake
#    1855 RbmC_standard
#      85 RbmC_truncated
#       5 RbmC_with_b-helix
#      26 RbmC_with_partial_M1M2
#      31 RbmC_w/o_M1M2

### 7) visualize gene tree
cd itol_tree/

# struct_type_annot.txt, strip
cp ../struct_type_template.txt struct_type_annot.txt
awk -F"\t" 'NR==FNR{a[$1]=$2;next}{print $1","a[$2]","$2} ' ../struct_type2col.mapping <(grep -wf ../all.ids ../gene_types.2.tsv) >> struct_type_annot.txt

zip -r RbmC_Bap1_gene.tree.zip gene.tree taxon_annot.txt struct_type_annot.txt
upload_tree RbmC_Bap1_gene.tree.zip # 997 leaves
# itol link: https://itol.embl.de/tree/128231272221561717447013
# root at midpoint; add color ranges
cd ..

### 8) relabel RbmC truncated genes and build the tree
# create relabel_RbmC.ids (RbmC clade 1 + RbmC w/o M1M2 in Vc) and exclude.ids (RbmC w/o M1M2 in Vc)
cd itol_tree/
grep -vf exclude.ids relabel_RbmC.ids| awk -F"\t" 'NR==FNR{a[$1];next}$1 in a{if($2!="RbmC_standard"){print $1"\tRbmC_truncated"}else{print $1"\t"$2}}!($1 in a){print $1"\t"$2}' - ../gene_types.2.tsv > ../gene_types.2.updated.tsv
cd ..

awk -F"\t" 'NR==FNR{a[$1]=$2}NR!=FNR{if($1 in a){print $2"\t"a[$1]}}' gene_types.2.updated.tsv repr2dup.mapping > gene_types.rest.updated.tsv
cat gene_types.2.updated.tsv gene_types.rest.updated.tsv > gene_types.all.updated.tsv # 4065
grep Bap1 gene_types.all.updated.tsv|cut -f1 |parallel -k ../scripts/get_faa_seq.sh {}  > Bap1_genes_2062_latest.fasta
# overwrite mapping file
paste <(cut -f1 gene_types.all.updated.tsv|cut -f1-2 -d_) gene_types.all.updated.tsv > gid_gene_type.RbmC_Bap1.mapping
python ../scripts/pickle_file.py gid_gene_type.RbmC_Bap1.mapping

cut -f3 gid_gene_type.RbmC_Bap1.mapping|sort |uniq -c
#    1784 Bap1_standard
#      29 Bap1_truncated
#      11 Bap1_wo_57aa_truncated
#     238 Bap1_w/o_loop
#       1 RbmC_fake
#    1855 RbmC_standard
#     103 RbmC_truncated
#       5 RbmC_with_b-helix
#      19 RbmC_with_partial_M1M2
#      20 RbmC_w/o_M1M2

cd itol_tree/
cp ../struct_type_template.txt struct_type_annot.2.txt
awk -F"\t" 'NR==FNR{a[$1]=$2;next}{print $1","a[$2]","$2} ' ../struct_type2col.mapping <(grep -wf ../all.ids ../gene_types.2.updated.tsv) >> struct_type_annot.2.txt

zip -r RbmC_Bap1_gene.tree.2.zip gene.tree taxon_annot.txt struct_type_annot.2.txt
upload_tree RbmC_Bap1_gene.tree.2.zip # 997 leaves
# itol link: https://itol.embl.de/tree/128231273328481717602697
cd ..

### 8) get Bap1 and loop-less Bap1 signal peptides
mkdir -p bap1_sig_pep
pep_fas=/gpfs/gsfs12/users/yangy34/projects/vibrio_YJ_v4/align_Bap1_RbmC/signal_peptide/singalP6.0_res/signal_peptide_aln.fasta
seqkit fx2tab $pep_fas | grep -wf <(grep Bap1_w/o_loop ../gene_types.2.tsv|cut -f1) - | awk -F"\t" '{print $1"\t"substr($2,15,17)substr($2,36,100)}' |sed "s/-//g" | seqkit tab2fx -w 0 > loop-less_Bap1_sig_pep.fasta
seqkit fx2tab $pep_fas | grep -wf <(grep Bap1_standard ../gene_types.all.tsv|cut -f1) - | awk -F"\t" '{print $1"\t"substr($2,15,100)}' | sed "s/-//g" | seqkit tab2fx -w 0 > std_Bap1_sig_pep.fasta

# Output:
# gene_types.all.updated.tsv
# gid_gene_type.RbmC_Bap1.mapping
