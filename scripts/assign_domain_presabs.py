import os, sys
from Bio import SeqIO

fa_path = sys.argv[1]

domain_intervals = {"beta-gamma-crystallin": (48, 243), "beta propeller (blades 1-5)": (244, 552), 
                    "beta prism C1_L": (553, 659), "57_aa":(660, 716), 
                    "beta prism C1_R":(717, 765), "beta-propeller (blades 6-8)":(766, 968), 
                    "beta prism C2":(969, 1160)}

print("\t".join(["Gene name"] + list(domain_intervals.keys())))
# given aligned sequences in fasta file, printout the absence and presence of domains
for rc in SeqIO.parse(fa_path, format="fasta"):
    name = rc.id
    seq = str(rc.seq)
    domain_presabs = []
    for dm in domain_intervals:
        start, end = domain_intervals[dm]
        each_dm = seq[start-1:end]
#        print(dm)
        nongap_perc = (len(each_dm) - each_dm.count("-"))/len(each_dm)*100
        if nongap_perc > 70:
            domain_presabs.append("1")
        elif 30 < nongap_perc <= 70:
            domain_presabs.append("0.5")
        else:
            domain_presabs.append("0")
    print(name+"\t"+"\t".join(domain_presabs))
