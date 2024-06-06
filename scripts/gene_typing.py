import os, sys
from Bio import SeqIO

fa_path = sys.argv[1] # the alignment curated/checked in Geneious
gene_lst_file = sys.argv[2] # the gene list of interest (do not have to be all present in fa_path)
gene_type_1_file = sys.argv[3] # initial typing result
output = sys.argv[4]

### Read gene list
gene_list = []
with open(gene_lst_file) as f:
    for line in f:
        gene_list.append(line.strip())

### Read typing file
gene_types_1 = {}
with open(gene_type_1_file) as f:
    for line in f:
        gene, type_1 = line.strip().split("\t")
        gene_types_1[gene] = type_1

beta_prop_beta_helix_lst = ["GCA_003544875.1_03596","GCA_002608565.1_00214","GCA_019670485.1_03306","GCA_019670025.1_03371","GCA_905175395.2_03458"]
RbmC_fake = ["GCA_002933855.1_04065"]

domain_intervals = {"beta-gamma-crystallin": (48, 243), "beta propeller (blades 1-5)": (244, 552), 
                    "beta prism C1_L": (553, 659), "57_aa":(660, 716), 
                    "beta prism C1_R":(717, 765), "beta-propeller (blades 6-8)":(766, 968), 
                    "beta prism C2":(969, 1160)}

# Name    Minimum Maximum Length
# RbmC_sec_sig    1       47      47
# beta_gamma_crystallin   48      243     196
# Bap1_sec_sig    218     243     26
# beta_propeller_blades_1-5       244     552     309
# beta_prism_B_C1_left    553     659     107
# 57aa_loop       660     716     57
# beta_prism_B_C1_right   717     765     49
# beta_propeller_blades_6-8       766     968     203
# beta_prism_C2   969     1160    192

# check
# l = []
# with open("d") as f:
#     for line in f:
#         l.append(line.strip())

### Store aligned sequences to dict
seq_dict = {}
for rc in SeqIO.parse(fa_path, format="fasta"):
    name = rc.id
    seq = str(rc.seq)
    seq_dict[name] = seq

### Get gene types ###
gene_types = {}
for gene in gene_list:
    gene_type = ""
    
    if gene in beta_prop_beta_helix_lst:
        gene_type = "RbmC_with_b-helix"
        gene_types[gene] = gene_type
        # print(gene+"\t"+gene_type)
        continue

    if gene in RbmC_fake:
        gene_type = "RbmC_fake"
        gene_types[gene] = gene_type
        # print(gene+"\t"+gene_type)
        continue
    
    domain_presabs = []
    seq = seq_dict[gene]
    
    # check
    # if gene in l:
    #     print(gene, seq.count("-"), len(seq), seq.count("-")/len(seq))
    
    if seq.count("-") > len(seq)*0.49:
        gene_types[gene] = gene_types_1[gene] + "_truncated"
        # print(gene+"\t"+gene_type)
        continue

    for dm in domain_intervals:
        start, end = domain_intervals[dm]
        each_dm = seq[start-1:end]
        nongap_perc = (len(each_dm) - each_dm.count("-"))/len(each_dm)*100
        if nongap_perc > 70:
            domain_presabs.append(1)
        elif 30 < nongap_perc <= 70:
            domain_presabs.append(0.5)
        else:
            domain_presabs.append(0)

    if domain_presabs[6]!=0:
        # beta_prism_C2 != 0: RbmC-like genes
        if domain_presabs[0] == 0:
            # M1M2 == 0
            gene_type = "RbmC_w/o_M1M2"
        elif domain_presabs[0] == 0.5:
            # M1M2 == 0.5
            gene_type = "RbmC_with_partial_M1M2"
        else:
            # M1M2 == 1
            gene_type = "RbmC_standard"
    else:
        if domain_presabs[0]!=0:
            # although beta_prism_2 is absent, but M1M2 == 1, this gene must be a RbmC-like gene
            if domain_presabs[0] == 0:
                # M1M2 == 0
                gene_type = "RbmC_w/o_M1M2"
            elif domain_presabs[0] == 0.5:
                # M1M2 == 0.5
                gene_type = "RbmC_with_partial_M1M2"
            else:
                gene_type = "RbmC_standard"

        # beta_prism_C2 == 0: Bap1-like genes
        elif domain_presabs[3] !=0:
             gene_type = "Bap1_standard"
        else:
            gene_type = "Bap1_w/o_loop"
    
    gene_types[gene] = gene_type
    # print(name, domain_presabs)
    # print(gene+"\t"+gene_type)

g = open(output, "w")
for gene in gene_types:
    g.write(gene+"\t"+gene_types[gene]+"\n")
g.close()
