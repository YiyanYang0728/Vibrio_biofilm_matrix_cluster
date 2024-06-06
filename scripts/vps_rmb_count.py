from BCBio import GFF
from collections import defaultdict
import sys

gid = sys.argv[1]
in_file = sys.argv[2] #"/data/yangy34/projects/vibrio_YJ_v5/VPS_cluster/genomes/run_PFF_KO_all_v4/repr_pff_res_checked/pffres.*.annot.gff"
vps_gene_type = {"vps-I": ["vpsA", "vpsB", "vpsD", "vpsE", "vpsF", "vpsI", "vpsJ", "vpsK"], "vps-II": ["vpsL", "vpsM", "vpsN", "vpsO"]}

in_handle = open(in_file)
rec_list = []
vps_cl = defaultdict(list)
rbm_genes = []
for rec in GFF.parse(in_handle):
    for feature in rec.features:
        gene = feature.qualifiers["Name"][0]
        if feature.qualifiers["ClusterID"][0] != "Cl_NA":
            cl_key = rec.id + "|" + feature.qualifiers["ClusterID"][0]
            if gene.startswith("vps"):
                vps_gene = gene
                vps_cl[cl_key].append(vps_gene)
        if "rbm" in gene.lower() and "_like" not in gene.lower() or "bap1" in gene.lower():
            rbm_genes.append(gene)
in_handle.close()

vps_1_genes = []
vps_2_genes = []
for cl_key in vps_cl:
    vps_1 = set(vps_cl[cl_key]).intersection(vps_gene_type["vps-I"])
    vps_2 = set(vps_cl[cl_key]).intersection(vps_gene_type["vps-II"])
    vps_1_genes.extend(vps_1)
    vps_2_genes.extend(vps_2)
vps_1_gene_ct = str(len(set(vps_1_genes)))
vps_2_gene_ct = str(len(set(vps_2_genes)))
print(gid+"\t"+vps_1_gene_ct+"\t"+vps_2_gene_ct+"\t"+",".join(set(rbm_genes)))