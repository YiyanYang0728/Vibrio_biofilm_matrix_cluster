### Step 8. Find prophage genomes
# Input:
cp ../Step3.annotate_RbmB_RbmA/itol_tree/leaves.list .
cut -f1-2 -d_ leaves.list |sort -u > gid.list
mkdir -p fna
cat gid.list | parallel cp ../Step1.annotate_cluster/data/{}.fna fna/

# Process:
# run VirSorter2
mkdir -p vs2_out
for gid in `cat gid.list`
do
    virsorter run -w vs2_out/${gid} -i fna/${gid}.fna --min-length 1000 -j 32 all
done

# run it for half a day
ls vs2_out/GCA_*/final-viral-boundary.tsv | parallel -k "cut -f1,4,5 {}|tail -n+2" > vs2_prophage_regions.tsv

cat leaves.list | parallel -k ../scripts/check_gene_in_prophage.sh {} > leaves_phage_asso.tsv

# Output:
# leaves_phage_asso.tsv