import sys
sys.path.append("../..")

from annotateOffs import *

seqs = {}
for row in iterTsvRows(open("ren2014.tab")):
    seqs[row.seq] = row.seq

newSeqs = extendSeqs(seqs, "dm3", 0, 3)

ofh = open("../ren2015.tab", "w")
ofh.write("origSeq\tseq\tmodFreq\tdm3Pos\n")
for row in iterTsvRows(open("ren2014.tab")):
    hits = newSeqs[row.seq]
    for newSeq, pos in hits:
        row = [row.seq, newSeq, row.modFreq, pos]
        ofh.write("\t".join(row)+"\n")

print "output written to %s" % ofh.name
