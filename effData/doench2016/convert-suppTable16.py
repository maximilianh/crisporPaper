from collections import defaultdict
import operator

from numpy import mean, log2, sum, array
import glob

goodGenes = set(["CCDC101", "CUL3", "HPRT1", "MED12", "NF1", "NF2", "TADA1", "TADA2B"])

fname = "suppTable16.txt"
print "infile %s" % fname
outFname = "../doench2016_hg19.guides.tab"
outFname2 = "../doench20166tg_hg19.guides.tab"
outFname3 = "../doench2016plx_hg19.guides.tab"
print "outfiles: %s, %s, %s" % (outFname, outFname2, outFname3)
ofh = open(outFname, "w")
ofh2 = open(outFname2, "w")
ofh3 = open(outFname3, "w")
ofh.write("guide\tseq\tmodFreq\n")
ofh2.write("guide\tseq\tmodFreq\n")
ofh3.write("guide\tseq\tmodFreq\n")
guideId = 0

for line in open(fname):
    if line.lower().startswith("construct") or line[0]=="\t":
        continue
    fields =  line.rstrip("\n").split("\t")
    gene = fields[4]
    if gene not in goodGenes:
        continue
    azdLfc = fields[34]
    sixTgLfc = fields[35]
    plxLfc = fields[36]
    seq = fields[0]

    name = "%s-%d" % (gene, guideId)
    guideId += 1

    row = [name, seq, str(azdLfc)]
    ofh.write("\t".join(row)+"\n")

    row = [name, seq, str(sixTgLfc)]
    ofh2.write("\t".join(row)+"\n")

    row = [name, seq, str(plxLfc)]
    ofh3.write("\t".join(row)+"\n")

