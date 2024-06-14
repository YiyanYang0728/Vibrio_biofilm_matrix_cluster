### Step 1. Annotate biofilm matrix cluster
# Input:
wd=<your path storing Prokka genome annotation and Kofam annotation output>
mkdir -p data/
# make sure all .faa, .fna, .gff. .kofam files are together in data/ folder

# ls ${wd}/faa/*|parallel realpath {}|parallel ln -s {} data/
# ls ${wd}/fna/*|parallel realpath {}|parallel ln -s {} data/
# ls ${wd}/gff/*|parallel realpath {}|parallel ln -s {} data/
# ls ${wd}/kofam/*|parallel realpath {}|parallel ln -s {} data/

# Process:
### Annotate with prokfunfind
mkdir -p run_PFF
cd run_PFF
ls data/*.kofam.tsv |parallel basename {} .kofam.tsv > ko_gid.list
awk -v pwd=$PWD '{print $0"\t"pwd"/data"}' ko_gid.list > genome_list.tsv
prokfunfind -f /data/yangy34/softwares/ProkFunFind/data/VcBMC_ko/config.yaml -p 10 -g genome_list.tsv -o repr_pff_res/pffres

# Output:
# repr_pff_res/pffres/pffres.*.annot.gff
