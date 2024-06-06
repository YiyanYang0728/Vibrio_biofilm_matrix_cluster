from Bio import SeqIO
import sys
records = SeqIO.parse(sys.argv[1], "fasta")
count = SeqIO.write(records, sys.argv[2], "clustal")
print("Converted %i records" % count)