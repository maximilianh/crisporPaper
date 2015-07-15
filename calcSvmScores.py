# collect various guide sequences and write their SVM scores to svmScores.tab
# so we don't have to calculate them each time, the R startup is rather slow

from annotateOffs import *
import os

seqs = set()

for fname in glob.glob("effData/*.tab"):
    print "reading %s" % fname
    for row in iterTsvRows(fname):
        seqs.add(row.seq[:20].upper())

scoreDict = calcSvmEffScores(seqs)
writeDict(scoreDict, "svmScores.tab")
logging.info("Wrote %d svm scores to %s" % (len(seqs), fname))
