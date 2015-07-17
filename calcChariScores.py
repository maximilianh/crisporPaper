# collect various guide sequences and write their SVM scores to chiariScores.tab
# so we don't have to calculate them each time, the SVM startup is rather slow

from annotateOffs import *
import os

seqs = set()

for fname in glob.glob("effData/*.tab"):
    print "reading %s" % fname
    for row in iterTsvRows(fname):
        assert(len(row.seq)==23)
        seqs.add(row.seq.upper())

scoreDict = calcChiariScores(seqs)
writeDict(scoreDict, "chiariScores2.tab")
logging.info("Wrote %d chiari scores to %s" % (len(seqs), fname))
