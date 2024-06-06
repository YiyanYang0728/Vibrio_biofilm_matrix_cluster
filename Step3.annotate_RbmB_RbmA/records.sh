### Step 3. Annotate RbmB and RbmA proteins

##### 3.1 annotate RbmB proteins
# get curated RbmC and Bap1 ProkFunFind results
mkdir -p repr_pff_res_RbmC_checked
scr=../scripts/curate_Bap1_RbmABC_VpsEF.py
cut -f1 -d. ../Step1.annotate_cluster/ko_gid.list > gid.list
cat gid.list | parallel -k python ${scr} {} ../Step1.annotate_cluster/repr_pff_res/*{}*.gff ../Step2.annotate_RbmC_Bap1/gid_gene_type.RbmC_Bap1.mapping.pkl ../Step1.annotate_cluster/data/{}.*.gff repr_pff_res_RbmC_checked/pffres.{}.annot.gff

# get gene predicted to have G3DSA:2.160.20.10 domain
cat gid.list | parallel ../scripts/get_rhbh_domain.sh {} | sort -u > rhbh_domain.list

# get putative RbmB ids
cat gid.list | parallel -k "grep 'Name=rbmB' repr_pff_res_RbmC_checked/pffres.{}.annot.gff|grep -oP 'ID=\KGCA_[0-9]*\.[0-9]*_[0-9]*'" > rbmB_ids.list

# get putative RbmB with G3DSA:2.160.20.10 domain
sort <(cut -f1 rhbh_domain.list|sort -u) <(cut -f1 rbmB_ids.list|sort -u)|uniq -d > rbmB_rhbh.list

# get gene name, length and if they are putative RbmB
cut -f1 rhbh_domain.list | parallel -k ../scripts/get_gene_name.sh {} > rhbh_domain_name.tsv
paste <(sort rhbh_domain.list|cut -f1,2) <(sort rhbh_domain_name.tsv|cut -f2) > rhbh_domain_name_len.tsv
awk -F"\t" 'NR==FNR{a[$1];next}$1 in a{print $1"\t"$2"\t"$3"\tRbmB"}!($1 in a){print $1"\t"$2"\t"$3"\tnon-RbmB"}' rbmB_rhbh.list rhbh_domain_name_len.tsv > rhbh_domain_name_len_label.tsv

# get more refined RbmB: a putative RbmB should be near either an RbmA or an RbmC; should be annotated as either rbmB hits based on hmm search or with G3DSA:2.160.20.10 domain
cat gid.list | parallel -k python ../scripts/get_rbmB_ids_putative.py {} repr_pff_res_RbmC_checked/pffres.{}.annot.gff ../Step1.annotate_cluster/data/{}.*.gff rhbh_domain.list > rbmB_ids_putative.tsv # 1750

# find split RbmB genes:
cut -f1 -d. rbmB_ids_putative.tsv|sort |uniq -d|grep -wf - rbmB_ids_putative.tsv|awk 'BEGIN{OFS="\t"}{split($1,a,".");print a[1],$1,$2}'|datamash -s -g 1 collapse 2 sum 3
# GCA_000175995   GCA_000175995.1_01868,GCA_000175995.1_01869     368
# GCA_000195065   GCA_000195065.1_00885,GCA_000195065.1_00886     416
# GCA_000223095   GCA_000223095.2_01747,GCA_000223095.2_01748     343
# GCA_000304915   GCA_000304915.2_01109,GCA_000304915.2_01110     303
# GCA_000696385   GCA_000696385.1_00045,GCA_000696385.1_00046     401
# GCA_000740515   GCA_000740515.2_02729,GCA_000740515.2_02730     427
# GCA_000966395   GCA_000966395.1_03744,GCA_000966395.1_03745     290
# GCA_001253315   GCA_001253315.1_01175,GCA_001253315.1_01176     398
# GCA_001402175   GCA_001402175.1_03452,GCA_001402175.1_03453     322
# GCA_001515165   GCA_001515165.1_02344,GCA_001515165.1_02345     285
# GCA_002073995   GCA_002073995.1_04201,GCA_002073995.1_04202     406
# GCA_003594515   GCA_003594515.1_01190,GCA_003594515.1_02369     307 # split by contigs
# GCA_003594945   GCA_003594945.1_01147,GCA_003594945.1_01148     354
# GCA_003716375   GCA_003716375.1_02732,GCA_003716375.1_02733     343
# GCA_008464965   GCA_008464965.1_01503,GCA_008464965.1_01504     403
# GCA_009763985   GCA_009763985.1_02732,GCA_009763985.1_02733     343
# GCA_009764005   GCA_009764005.1_00989,GCA_009764005.1_00990     402
# GCA_013357785   GCA_013357785.1_01930,GCA_013357785.1_01931     354

cut -f1 -d. rbmB_ids_putative.tsv|sort|uniq -u|grep -wf - rbmB_ids_putative.tsv | awk 'BEGIN{OFS="\t"}{split($1,a,".");print a[1],$1,$2}' > rbmB_gid2ids_putative.tsv
cut -f1 -d. rbmB_ids_putative.tsv|sort |uniq -d|grep -wf - rbmB_ids_putative.tsv|awk 'BEGIN{OFS="\t"}{split($1,a,".");print a[1],$1,$2}'|datamash -s -g 1 collapse 2 sum 3 >> rbmB_gid2ids_putative.tsv
# 1732

# add putative genes to rhbh_domain_name_len_label.tsv ==> rhbh_domain_name_len_label_2.tsv
awk -F"\t" 'NR==FNR{a[$1];next}$1 in a{print $1"\t"$2"\t"$3"\tRbmB"}!($1 in a){print $1"\t"$2"\t"$3"\tnon-RbmB"}' <(cut -f2 rbmB_gid2ids_putative.tsv|tr "," "\n") rhbh_domain_name_len.tsv > rhbh_domain_name_len_label_2.tsv
sort <(cut -f2 rbmB_gid2ids_putative.tsv|tr "," "\n") <(cut -f1 rhbh_domain_name_len.tsv) <(cut -f1 rhbh_domain_name_len.tsv)|uniq -u |parallel ../scripts/get_gene_len.sh {} > diff_ids.len
awk '{print $1"\t"$2"\thypothetical protein\tRbmB"}' diff_ids.len >> rhbh_domain_name_len_label_2.tsv

# add seq similar genes to rhbh_domain_name_len_label_2.tsv ==> rhbh_domain_name_len_label_3.tsv
awk '$NF=="non-RbmB"{print $1}' rhbh_domain_name_len_label_2.tsv |parallel -k ../scripts/get_faa_seq.sh {} > non_RbmB_genes.fasta
awk '$NF=="RbmB"{print $1}' rhbh_domain_name_len_label_2.tsv |parallel -k ../scripts/get_faa_seq.sh {} > RbmB_genes.fasta
makeblastdb -in non_RbmB_genes.fasta -dbtype prot
blastp -query RbmB_genes.fasta -db non_RbmB_genes.fasta -num_threads 32 -outfmt 6 > RbmB_vs_nonRbmB_blastp.out
awk -F"\t" 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;next}{print $0"\t"a[$1]"\t"a[$2]}' rhbh_domain_name_len_label_2.tsv RbmB_vs_nonRbmB_blastp.out > RbmB_vs_nonRbmB_blastp_qrylen_sbjlen.out
awk -F"\t" '$3>60 && ($4/$NF>0.9 || $4/$(NF-1)>0.9)' RbmB_vs_nonRbmB_blastp_qrylen_sbjlen.out|cut -f2|sort -u|awk -F"\t" 'BEGIN{OFS="\t"}NR==FNR{a[$1];next}$1 in a{print $1,$2,$3,"RbmB"}!($1 in a){print }' - rhbh_domain_name_len_label_2.tsv > rhbh_domain_name_len_label_3.tsv

awk '$NF=="RbmB"' rhbh_domain_name_len_label_2.tsv |cut -f1 > rhbh_domain_name_len_label_2.rbmB.ids
awk '$NF=="RbmB"' rhbh_domain_name_len_label_3.tsv |cut -f1 > rhbh_domain_name_len_label_3.rbmB.ids
sort rhbh_domain_name_len_label_2.rbmB.ids rhbh_domain_name_len_label_2.rbmB.ids rhbh_domain_name_len_label_3.rbmB.ids | uniq -u > extra_rbmB_genes
# # 10 extra genes are added according to seq similarity:
# GCA_000287115.2_01678
# GCA_000287155.2_00286
# GCA_001252875.1_00016
# GCA_001253155.1_00992
# GCA_001860345.1_03538
# GCA_001860485.1_01620
# GCA_002204085.1_03417
# GCA_002204105.1_03460
# GCA_016456785.1_01465
# GCA_017426745.1_02752

cut -f1 -d. extra_rbmB_genes |grep -wf - ../Step2.annotate_RbmC_Bap1/gid_gene_type.RbmC_Bap1.mapping|cut -f1|sort -u|wc -l
# 10
# all these genes have RbmC/Bap1 in their genomes, so they are likely the true RbmB genes and are saparate from RbmC due to contig fragmentation

# create mapping file
awk '{split($1,a,"_"); print a[1]"_"a[2]"\t"$1"\tRbmB"}' rhbh_domain_name_len_label_3.rbmB.ids > rbmB.mapping
cat rbmB.mapping ../Step2.annotate_RbmC_Bap1/gid_gene_type.RbmC_Bap1.mapping > gid_gene_type.RbmBC_Bap1.mapping
python ../scripts/pickle_file.py gid_gene_type.RbmBC_Bap1.mapping

awk '$2<1000 && $2>=200{print $1}' rhbh_domain_name_len_label_3.tsv | parallel -k ../scripts/get_faa_seq.sh {} > rhbh_domain_len_filtered.fasta
seqkit rmdup -s rhbh_domain_len_filtered.fasta -D rhbh_domain_len_filtered.dup.list -d rhbh_domain_len_filtered.dup.tsv -o rhbh_domain_len_filtered_rmdup.fasta
# 2545 unique seqs left

# align seqs seq_len ==> 4584 (too long)
mafft --dash --thread 32 --localpair rhbh_domain_len_filtered_rmdup.fasta > rhbh_domain_len_filtered_rmdup.mafft.aln
trimal -gt 0.2 -in rhbh_domain_len_filtered_rmdup.mafft.aln -out trimal_out.fasta
# MSA seq_len ==> 1017
seqkit fx2tab trimal_out.fasta | awk '$1~/GCA/' | seqkit tab2fx -w 0 > rhbh_domain_rmdup_trimmed.fasta

##### 3.2 annotate RbmA proteins
# get curated RbmC, Bap1 and RbmB ProkFunFind results
mkdir -p repr_pff_res_RbmBC_checked
scr=../scripts/curate_Bap1_RbmABC_VpsEF.py
cat gid.list | parallel -k python ${scr} {} repr_pff_res_RbmC_checked/pffres.{}.annot.gff gid_gene_type.RbmBC_Bap1.mapping.pkl ../Step1.annotate_cluster/data/{}.*.gff repr_pff_res_RbmBC_checked/pffres.{}.annot.gff

# get gene predicted to have  domain
cat gid.list | parallel ../scripts/get_FnIII_domain.sh {} > FnIII_domain.list

# get putative RbmA ids
cat gid.list | parallel -k "grep 'Name=rbmA' repr_pff_res_RbmBC_checked/pffres.{}.annot.gff|grep -oP 'ID=\KGCA_[0-9]*\.[0-9]*_[0-9]*'" | sort -u > rbmA_ids.list

# get putative RbmB with G3DSA:2.60.40.3880 domain
sort <(cut -f1 FnIII_domain.list|sort -u) <(cut -f1 rbmA_ids.list|sort -u)|uniq -d > rbmA_FnIII.list

# get gene name, length and if they are putative RbmA
cut -f1 FnIII_domain.list | parallel -k ../scripts/get_gene_name.sh {} > FnIII_domain_name.tsv
paste <(sort FnIII_domain.list|cut -f1,2) <(sort FnIII_domain_name.tsv|cut -f2) | sort -u > FnIII_domain_name_len.tsv
awk -F"\t" 'NR==FNR{a[$1];next}$1 in a{print $1"\t"$2"\t"$3"\tRbmA"}!($1 in a){print $1"\t"$2"\t"$3"\tnon-RbmA"}' rbmA_FnIII.list FnIII_domain_name_len.tsv > FnIII_domain_name_len_label.tsv

# get more refined RbmA: a putative RbmA should be near either an RbmB or an RbmC; should be annotated as either rbmA hits based on hmm search or with G3DSA:2.60.40.3880 domain
cat gid.list | parallel -k python ../scripts/get_rbmA_ids_putative.py {} repr_pff_res_RbmBC_checked/pffres.{}.annot.gff ../Step1.annotate_cluster/data/{}.*.gff FnIII_domain.list > rbmA_ids_putative.tsv

# find split RbmA genes:
cut -f1 -d. rbmA_ids_putative.tsv|sort |uniq -d|grep -wf - rbmA_ids_putative.tsv|awk 'BEGIN{OFS="\t"}{split($1,a,".");print a[1],$1,$2}'|datamash -s -g 1 collapse 2 sum 3
cut -f1 -d. rbmA_ids_putative.tsv|sort |uniq -d > rbmA_gids_split.list
# 13 RbmA genes are split
# GCA_000237745   GCA_000237745.2_01765,GCA_000237745.2_01766     258
# GCA_000754625   GCA_000754625.1_02797,GCA_000754625.1_02798     265
# GCA_000786345   GCA_000786345.1_02787,GCA_000786345.1_02788     268
# GCA_000939665   GCA_000939665.1_02681,GCA_000939665.1_02682     265
# GCA_000966375   GCA_000966375.1_02604,GCA_000966375.1_02605     268
# GCA_000966385   GCA_000966385.1_02943,GCA_000966385.1_02944     265
# GCA_001402175   GCA_001402175.1_03454,GCA_001402175.1_03455,GCA_001402175.1_03456       268
# GCA_001641745   GCA_001641745.1_02237,GCA_001641745.1_02238     265
# GCA_001857145   GCA_001857145.1_02312,GCA_001857145.1_02313     258
# GCA_001857245   GCA_001857245.1_02333,GCA_001857245.1_02334     258
# GCA_009762985   GCA_009762985.1_01254,GCA_009762985.1_01255     258
# GCA_009763045   GCA_009763045.1_00367,GCA_009763045.1_00368     255
# GCA_009763105   GCA_009763105.1_02894,GCA_009763105.1_02895     263

# check RbmA genes only having one FnIII domain
grep -vf rbmA_gids_split.list rbmA_ids_putative.tsv | cut -f1|grep -wf - FnIII_domain.list|cut -f1|sort |uniq -c|awk '$1==1{print $2}' > rbmA_ids_only_1_domain.list
grep -wf rbmA_ids_only_1_domain.list rbmA_ids_putative.tsv
# too short to have 2 FnIII domains
# GCA_001250615.1_03069   69.0    rbmB    GCA_001250615.1_03070   Cl_0    1       rbmC    GCA_001250615.1_03071   Cl_0    2       fniii   rbma
# GCA_008083545.1_01663   87.0    rbmB    GCA_008083545.1_01662   Cl_0    -1      rbmC    GCA_008083545.1_01661   Cl_0    -2      fniii   rbma
# GCA_013155105.1_03686   166.0   rbmB    GCA_013155105.1_03687   Cl_0    1       rbmC    -       -       -       fniii   rbma

grep -vf <(cut -f1 rbmA_ids_putative.tsv) FnIII_domain.list | cut -f1 |sort |uniq -c|awk '$1>=2 {print $2}' > 151_nonRbmA_2_FnIII_ids.list
# 151 genes are not rbmA genes but have >=2 FnIII domains and all of these genes only have 2 FnIII domains; check their taxonomy and RbmC presence, they are all likely to be real RbmA

cat <(cut -f1 rbmA_ids_putative.tsv) 151_nonRbmA_2_FnIII_ids.list > rbmA_ids_putative_2.list

# check other FnIII domain genes
grep -vf rbmA_ids_putative_2.list FnIII_domain.list | cut -f1 > 68_nonRbmA_1_FnIII_ids.list

grep -wf 68_nonRbmA_1_FnIII_ids.list rbmA_FnIII.list|cut -f1|sort -u
# check that 5 genes are rbmA hmm hits and with 1 FnIII domain:
# GCA_000473785.1_00868: RbmA, split with GCA_000473785.1_00869, 97, edge of contig
# GCA_001249315.1_03251: RbmA, split with GCA_001249315.1_03252, 37
# GCA_001249315.1_03252: RbmA, split with GCA_001249315.1_03251, 186
# GCA_001249795.1_03493: RbmA, not-split, 140
# GCA_015343015.1_02539: maybe RbmA, not-split, far from rbm cluster, 95, edge of contig

# create a file to store these genes: 5_nonRbmA_1_FnIII_ids.list and add to rbmA list
cat 5_nonRbmA_1_FnIII_ids.list rbmA_ids_putative_2.list > rbmA_ids_putative_3.list

# check the rest in 68 genes
grep -vf 5_nonRbmA_1_FnIII_ids.list 68_nonRbmA_1_FnIII_ids.list > 63_nonRbmA_1_FnIII_ids.list
# non of these ids have checked rbmC nearby; check their taxon that they are in s__Vibrio campbellii, s__Vibrio diabolicus and s__Vibrio parahaemolyticus; 
# by looking at context, they are **less likely** to be real RbmA genes.

# create mapping file
awk '{split($1,a,"_"); print a[1]"_"a[2]"\t"$1"\tRbmA"}' rbmA_ids_putative_3.list > rbmA.mapping
cat rbmA.mapping gid_gene_type.RbmBC_Bap1.mapping > gid_gene_type.RbmABC_Bap1.mapping
python ../scripts/pickle_file.py gid_gene_type.RbmABC_Bap1.mapping

# Output:
# gid_gene_type.RbmBC_Bap1.mapping
# gid_gene_type.RbmABC_Bap1.mapping
