# determine efficiency score cutoffs for the 75 percentile

from annotateOffs import *
from collections import defaultdict
import numpy as np

scores = defaultdict(list)
for fname in glob.glob("effData/*.scores.tab"):
    if "hart" in fname and "Avg" not in fname:
        continue
    print "reading", fname
    for row in iterTsvRows(fname):
        for scoreType, score in zip(row._fields[7:], row[7:22]):
            score = float(score)
            # make sure doench score are not the old ones in the range 0-1.0
            if scoreType=="doench" and score < 1.0 and score!=-1 and score!=0.0:
                print fname, score, scoreType
                assert(False)
            #if scoreType=="doench":
                #print fname, scoreType, score
            scores[scoreType].append(float(score))

print "got scores for %d sequences" % len(scores["doench"])

# transform to scoreType -> list of scores
#effDicts = calcEffScores(seqs)
#for seq, effDict in effDicts.iteritems():
    #for scoreType in getScoreTypes():
        #scores[scoreType].append(effDict[scoreType])

cutoffs = {}
for scoreType, scoreList in scores.items():
    scoreList.sort()
    cutoffs[scoreType] = np.percentile(scoreList, 75)
print cutoffs
