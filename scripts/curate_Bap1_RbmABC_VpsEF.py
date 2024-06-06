from BCBio import GFF
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
import sys, pickle

gid = sys.argv[1]
in_file = sys.argv[2] #"/data/yangy34/projects/vibrio_YJ_v5/VPS_cluster/genomes/run_PFF_KO_all_v4/repr_pff_res/pffres.GCA_006802795.1.annot.gff"
mapping_file = sys.argv[3] # gid_gene_type.VpsEF_RbmBC_Bap1.mapping.pkl
gff_file = sys.argv[4] #/data/yangy34/projects/vibrio_YJ_v5/VPS_cluster/genomes/gff/GCA_006802795.1.gff
outfile = sys.argv[5] #"/data/yangy34/test/GCA_006802795.gff"
lowercs = lambda s: s[:1].lower() + s[1:] if s else ''

# read dict
if not mapping_file.endswith(".pkl"):
    gene2type_all = defaultdict()
    with open(mapping_file) as f:
        for line in f:
            cols = line.strip().split("\t")
            gid, gene, typ = cols
            gene2type_all[gene] = typ
else:
    f = open(mapping_file, "rb")
    gene2type_all = pickle.load(f)

search_genes = list(set([lowercs(i.split("_")[0]) for i in list(set(gene2type_all.values()))]))
# print(search_genes)

# in_handle = open(in_file)
gff_features = defaultdict()
for rec in GFF.parse(gff_file):
    # get this genome's gene2type
    gene2type = {gene:gene2type_all[gene] for gene in gene2type_all if gene.split(".")[0]==rec.id.split(".")[0]}
    for feature in rec.features:
        if feature.id in gene2type:
            gff_features[feature.id] = (rec, feature)
# in_handle.close()

### what if there is any gene not in original pff res but are annotated as vpsEF and RbmBC?
in_handle = open(in_file)
rec_list = []
seen_genes = []
for rec in GFF.parse(in_handle):
    for feature in rec.features:
        seen_genes.append(feature.id)
        if feature.id not in gene2type:
            if feature.qualifiers['Name'][0] in search_genes:
                feature.qualifiers['Name'] = [feature.qualifiers['Name'][0]+"_like"]
            #     feature.qualifiers['Name'] = ['rbmC_like']
            # if feature.qualifiers['Name'] == ["rbmC"]:
            #     feature.qualifiers['Name'] = ['rbmC_like']
            # if feature.qualifiers['Name'] == ["vpsE"]:
            #     feature.qualifiers['Name'] = ['vpsE_like']
            # if feature.qualifiers['Name'] == ["vpsF"]:
            #     feature.qualifiers['Name'] = ['vpsF_like']
            # if feature.qualifiers['Name'] == ["rbmB"]:
            #     feature.qualifiers['Name'] = ['rbmB_like']
            # if feature.qualifiers['Name'] == ["rbmA"]:
            #     feature.qualifiers['Name'] = ['rbmA_like']
        else:
            feature.qualifiers['Name'] = [lowercs(gene2type[feature.id])] # update name
            feature.qualifiers['Status'] = ["checked"] # add Status=checked
    rec_list.append(rec)
in_handle.close()

for gene in (set(gff_features.keys()) - set(seen_genes)):
    rec, feature = gff_features[gene]
    # update feature
    feature.qualifiers['source'] = ['ProkFunFind']
    feature.qualifiers['ClusterID'] = ['Cl_NA']
    feature.qualifiers['Name'] = [gene2type[gene]]
    feature.qualifiers['Status'] = ["checked"]
    feature.qualifiers['Target'] = ["None"]
    feature.qualifiers['evalue'] = ["0.0"]
    feature.qualifiers.pop('inference', None)
    feature.qualifiers.pop('locus_tag', None)
    feature.qualifiers.pop('product', None)
    for record in rec_list:
        # add feature to existing rec
        if record.id == rec.id:
            record.features.append(feature)
            break
    else:
        # create a rec
        rec.annotations={}
        rec.features= [feature]
        rec_list.append(rec)

# clean rec.annotations
for rec in rec_list:
    rec.annotations = {}
  
with open(outfile, "w") as out_handle:
    GFF.write(rec_list, out_handle)
