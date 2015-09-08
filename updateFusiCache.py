# write all fusi scores to a cache file

import pickle
scores = {}
from annotateOffs import *

for fname in glob.glob("effData/*.scores.tab"):
    print "reading %s" % fname
    for row in iterTsvRows(fname):
        seq = row.longSeq100Bp
        score = float(row.fusi)
        if seq in scores:
            if scores[seq]!=score:
                assert((scores[seq]-score) < 0.01)
                #print seq, score
        #assert(seq not in scores or scores[seq]==score)
        scores[seq] = score

print "read %d scores" % len(scores)

writeDict(scores, "out/fusiCache.tab")
print "wrote to out/fusiCache.tab"
