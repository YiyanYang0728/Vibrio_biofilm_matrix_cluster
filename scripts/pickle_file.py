import sys, pickle, os

infile=sys.argv[1]
outfile=infile+".pkl"

d = {}
with open(infile) as f:
    for line in f:
        cols = line.strip().split("\t")
        d[cols[1]] = cols[2]

with open(outfile, "wb") as g:
    pickle.dump(d, g)
    