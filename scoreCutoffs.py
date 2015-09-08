# determine efficiency score cutoffs for the 75 percentile

from annotateOffs import *
from collections import defaultdict
import numpy as np

seqs = set()
for fname in glob.glob("effData/*.ext.tab"):
    #print fname
    for row in iterTsvRows(fname):
        seq = row.extSeq
        assert(len(seq)==34)
        seqs.add(seq)

print "got scores for %d sequences" % len(seqs)

# transform to scoreType -> list of scores
scores = defaultdict(list)
effDicts = calcEffScores(seqs, skipOof=True)
for seq, effDict in effDicts.iteritems():
    for scoreType in ["chariRank", "ssc", "svm", "doench"]:
        scores[scoreType].append(effDict[scoreType])

cutoffs = {}
for scoreType, scoreList in scores.items():
    cutoffs[scoreType] = np.percentile(scoreList, 75)
print cutoffs
