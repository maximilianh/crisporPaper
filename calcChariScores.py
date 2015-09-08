# collect various guide sequences and write their SVM scores to chariScores.tab
# so we don't have to calculate them each time, the SVM startup is rather slow

from annotateOffs import *
import os

seqs = set()

oldScores = readDict("chariScores.tab")

for fname in glob.glob("effData/*.tab"):
    print "reading %s" % fname
    for row in iterTsvRows(fname):
        if not len(row.seq)==23:
            print row
        assert(len(row.seq)==23)
        seq = row.seq.upper()
        if seq not in oldScores and "N" not in seq:
            seqs.add(seq)

print "%d sequences to do" % len(seqs)

for s in seqs:
    assert(len(s)==23)

scoreDict = calcChariScores(seqs)
scoreDict.update(oldScores)
writeDict(scoreDict, "chariScores.tab")
logging.info("Wrote %d chari scores to %s" % (len(seqs), fname))
