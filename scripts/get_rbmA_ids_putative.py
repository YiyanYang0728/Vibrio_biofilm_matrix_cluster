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
fniii_domain = sys.argv[4]
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

# read G3DSA:2.60.40.3880 list
fniii_id_list = []
with open(fniii_domain) as f:
    for line in f:
        cols = line.strip().split()
        if cols[0] not in fniii_id_list:
            fniii_id_list.append(cols[0])

in_handle = open(in_file)
rec_list = []
rbmB_list = []
rbmC_list = []
rbmA_id_list = []
vps_cl = defaultdict(list)
for rec in GFF.parse(in_handle):
    for feature in rec.features:
        if "ID" in feature.qualifiers:
            Id = feature.qualifiers["ID"][0]
            gene = feature.qualifiers["Name"][0]
            clstr = feature.qualifiers["ClusterID"][0]
            start = feature.location.start
            end = feature.location.end
            # only print out RbmA gene near a checked RbmC or checked RbmB ==> in a cluster
            if "Status" in feature.qualifiers:
                status = feature.qualifiers["Status"][0]
                if "rbmc" in gene.lower() and gene.lower()!="rbmc_fake":
                    rbmC = Gene(rec.id, clstr, Id, gene, start, end)
                    rbmC_list.append(rbmC)
                if "rbmb" in gene.lower():
                    rbmB = Gene(rec.id, clstr, Id, gene, start, end)
                    rbmB_list.append(rbmB)
            if "rbma" in gene.lower():
                rbmA = Gene(rec.id, clstr, Id, gene, start, end)
                rbmA_id_list.append(rbmA.Id)
in_handle.close()

picked_rbmA_dict = defaultdict(defaultdict)
other_gene = None
for q_gene in gene_list:
    for rbmC in rbmC_list:
        if rbmC!="" and q_gene.ctg == rbmC.ctg and q_gene.Id:
            other_gene = rbmC
            rbmC_dist = q_gene.gene_dist(rbmC)
            # print(rbmC_dist)
            # a putative RbmA is close to RbmC and have fniii domain
            if abs(rbmC_dist) <= 10 and q_gene.Id in fniii_id_list:
                picked_rbmA_dict[q_gene]["rbmC"] = (other_gene, rbmC_dist)
                break
            else:
                if abs(rbmC_dist) <= 10 and q_gene.Id in rbmA_id_list:
                    picked_rbmA_dict[q_gene]["rbmC"] = (other_gene, rbmC_dist)
                    break
    for rbmB in rbmB_list:
        if rbmB!="" and q_gene.ctg == rbmB.ctg and q_gene.Id in fniii_id_list:
            other_gene = rbmB
            rbmB_dist = q_gene.gene_dist(rbmB)
            # print(rbmB_dist)
            # a putative RbmA is close to RbmB and have fniii domain
            if abs(rbmB_dist) <= 10:
                picked_rbmA_dict[q_gene]["rbmB"] = (other_gene, rbmB_dist)
                break
            else:
                if abs(rbmB_dist) <= 10 and q_gene.Id in rbmA_id_list:
                    picked_rbmA_dict[q_gene]["rbmB"] = (other_gene, rbmB_dist)
                    break

# print([i.Id for i in picked_rbmA_dict])
for rbmA in picked_rbmA_dict:
    if "rbmB" in picked_rbmA_dict[rbmA]:
        rbmB, rbmB_dist = picked_rbmA_dict[rbmA]["rbmB"]
        rbmB_Id = rbmB.Id
        rbmB_clstr = rbmB.clstr
    else:
        rbmB_Id = "-"
        rbmB_clstr = "-"
        rbmB_dist = "-"
        
    if "rbmC" in picked_rbmA_dict[rbmA]:
        rbmC, rbmC_dist = picked_rbmA_dict[rbmA]["rbmC"]
        rbmC_Id = rbmC.Id
        rbmC_clstr = rbmC.clstr
    else:
        rbmC_Id = "-"
        rbmC_clstr = "-"
        rbmC_dist = "-"
        
    if rbmA.Id in fniii_id_list:
        label1 = "fniii"
    else:
        label1 = "-"
    if rbmA.Id in fniii_id_list:
        label2 = "rbma"
    else:
        label2 = "-"
    row = [rbmA.Id, str(rbmA.length), "rbmB", rbmB_Id, rbmB_clstr, str(rbmB_dist), "rbmC", rbmC_Id, rbmC_clstr, str(rbmC_dist), label1, label2]
    print("\t".join(row))
