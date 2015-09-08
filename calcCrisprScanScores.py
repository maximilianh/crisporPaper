# collect all guide sequences and write their crisprScan scores to a smaller table
# so we don't have to iterate over the huge file each time

from annotateOffs import *
import os

seqs = set()

oldScores = readDict("crisprScanScores.tab")

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

scoreDict = {}
doneSeqs = set()
for line in open("crisprScan/allScores.tab"):
    seq, score = line.rstrip("\n").split("\t")
    if seq in seqs:
        scoreDict[seq] = score
        doneSeqs.add(seq)

notDoneSeqs = (seqs-doneSeqs)
print "%d seqs not found out of %d" % (len(notDoneSeqs), len(seqs))
scoreDict.update(oldScores)
writeDict(scoreDict, "crisprScanScores.tab")
logging.info("Wrote %d crisprScan scores to %s" % (len(seqs), fname))
