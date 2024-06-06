from BCBio import GFF
from collections import defaultdict
import sys

class Gene:
    def __init__(self, ctg, clstr, Id, gene, start, end):
        self.ctg = ctg
        self.clstr = clstr
        self.Id = Id
        self.gene = gene
        self.length = (end-start)/3
    def gene_dist(self, other_gene):
        return(int(other_gene.Id[-5:])-int(self.Id[-5:]))


gid = sys.argv[1]
in_file = sys.argv[2] # "repr_pff_res_checked/pffres.*.annot.gff"
gff_file = sys.argv[3] # "*.gff"
rhbh_domain = sys.argv[4]
gene_list = []
in_handle = open(gff_file)
for rec in GFF.parse(in_handle):
    for feature in rec.features:
        if "ID" in feature.qualifiers:
            Id = feature.qualifiers["ID"][0]
            gene = feature.qualifiers["Name"][0] if "Name" in feature.qualifiers else None
            clstr = None
            start = feature.location.start
            end = feature.location.end
            query_gene = Gene(rec.id, clstr, Id, gene, start, end)
            gene_list.append(query_gene)
in_handle.close()

# read G3DSA:2.160.20.10 list

rhbh_id_list = []
with open(rhbh_domain) as f:
    for line in f:
        cols = line.strip().split()
        rhbh_id_list.append(cols[0])

in_handle = open(in_file)
rec_list = []
rbmA_list = []
rbmC_list = []
rbmB_id_list = []
vps_cl = defaultdict(list)
for rec in GFF.parse(in_handle):
    for feature in rec.features:
        if "ID" in feature.qualifiers:
            Id = feature.qualifiers["ID"][0]
            gene = feature.qualifiers["Name"][0]
            clstr = feature.qualifiers["ClusterID"][0]
            start = feature.location.start
            end = feature.location.end
            # only print out RbmB gene near a checked RbmC ==> in a cluster
            if "Status" in feature.qualifiers:
                status = feature.qualifiers["Status"][0]
                if "rbmc" in gene.lower() and gene.lower()!="rbmc_fake":
                    rbmC = Gene(rec.id, clstr, Id, gene, start, end)
                    rbmC_list.append(rbmC)
            if "rbma" in gene.lower():
                rbmA = Gene(rec.id, clstr, Id, gene, start, end)
                rbmA_list.append(rbmA)
            if "rbmb" in gene.lower():
                rbmB = Gene(rec.id, clstr, Id, gene, start, end)
                rbmB_id_list.append(rbmB.Id)
in_handle.close()
# print(rbmB_id_list)

picked_rbmB_dict = defaultdict(defaultdict)
other_gene = None
for q_gene in gene_list:
    for rbmC in rbmC_list:
        if rbmC!="" and q_gene.ctg == rbmC.ctg and q_gene.Id:
            other_gene = rbmC
            rbmC_dist = q_gene.gene_dist(rbmC)
            # print(rbmC_dist)
            # a putative RbmB is close to RbmC and have rhbh domain
            if abs(rbmC_dist) <= 8 and q_gene.Id in rhbh_id_list:
                picked_rbmB_dict[q_gene]["rbmC"] = (other_gene, rbmC_dist)
                break
            else:
                if abs(rbmC_dist) <= 8 and q_gene.Id in rbmB_id_list:
                    picked_rbmB_dict[q_gene]["rbmC"] = (other_gene, rbmC_dist)
                    break
    for rbmA in rbmA_list:
        if rbmA!="" and q_gene.ctg == rbmA.ctg and q_gene.Id in rhbh_id_list:
            other_gene = rbmA
            rbmA_dist = q_gene.gene_dist(rbmA)
            # print(rbmA_dist)
            # a putative RbmB is close to RbmA and have rhbh domain
            if abs(rbmA_dist) <= 8:
                picked_rbmB_dict[q_gene]["rbmA"] = (other_gene, rbmA_dist)
                break
            else:
                if abs(rbmA_dist) <= 8 and q_gene.Id in rbmB_id_list:
                    picked_rbmB_dict[q_gene]["rbmA"] = (other_gene, rbmA_dist)
                    break

# print([i.Id for i in picked_rbmB_dict])
for rbmB in picked_rbmB_dict:
    if "rbmA" in picked_rbmB_dict[rbmB]:
        rbmA, rbmA_dist = picked_rbmB_dict[rbmB]["rbmA"]
        rbmA_Id = rbmA.Id
        rbmA_clstr = rbmA.clstr
    else:
        rbmA_Id = "-"
        rbmA_clstr = "-"
        rbmA_dist = "-"
        
    if "rbmC" in picked_rbmB_dict[rbmB]:
        rbmC, rbmC_dist = picked_rbmB_dict[rbmB]["rbmC"]
        rbmC_Id = rbmC.Id
        rbmC_clstr = rbmC.clstr
    else:
        rbmC_Id = "-"
        rbmC_clstr = "-"
        rbmC_dist = "-"
        
    if rbmB.Id in rhbh_id_list:
        label1 = "rhbh"
    else:
        label1 = "-"
    if rbmB.Id in rbmB_id_list:
        label2 = "rbmb"
    else:
        label2 = "-"
    row = [rbmB.Id, str(rbmB.length), "rbmA", rbmA_Id, rbmA_clstr, str(rbmA_dist), "rbmC", rbmC_Id, rbmC_clstr, str(rbmC_dist), label1, label2]
    print("\t".join(row))
