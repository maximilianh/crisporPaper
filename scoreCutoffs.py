# determine efficiency score cutoffs for the 75 percentile

from annotateOffs import *
from collections import defaultdict
import numpy as np

scores = defaultdict(list)
for fname in glob.glob("effData/*.scores.tab"):
    print "reading", fname
    for row in iterTsvRows(fname):
        for scoreType, score in zip(row._fields[7:], row[7:]):
            scores[scoreType].append(float(score))

print "got scores for %d sequences" % len(scores["doench"])

# transform to scoreType -> list of scores
#effDicts = calcEffScores(seqs)
#for seq, effDict in effDicts.iteritems():
    #for scoreType in getScoreTypes():
        #scores[scoreType].append(effDict[scoreType])

cutoffs = {}
for scoreType, scoreList in scores.items():
    print scoreList
    cutoffs[scoreType] = np.percentile(scoreList, 75)
print cutoffs
