#create a table with all scores for all sequences in all datasets
import glob
from annotateOffs import *
from collections import defaultdict
from os.path import basename

scoreTypes = ["chariRank", "ssc", "svm", "doench"]

headers = None
outRows = []
for fname in glob.glob("effData/*.ext.tab"):
    dataset = basename(fname).split(".")[0]
    seqs = set()
    #print fname
    for row in iterTsvRows(fname):
        headers = row._fields
        seq = row.extSeq
        assert(len(seq)==34)
        seqs.add(seq)

    effDicts = calcEffScores(seqs, skipOof=True)
    print "scored for %d sequences" % len(seqs)

    for row in iterTsvRows(fname):
        newRow = list(row)
        newRow.insert(0, dataset)
        effDict = effDicts[row.extSeq]
        for scoreType in scoreTypes:
            newRow.append(effDict[scoreType])
        outRows.append(newRow)

headers = list(headers)
headers.insert(0, "dataset")
headers.extend(scoreTypes)

ofh = open("out/effScores.tab", "w")
ofh.write("\t".join(headers)+"\n")
for row in outRows:
    row = [str(x) for x in row]
    ofh.write("\t".join(row)+"\n")
print "wrote to %s" % ofh.name

