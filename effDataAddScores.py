# add all efficiency scores to the .context.tab files in effData/
# and write to .scores.tab files

from annotateOffs import *
import glob
import sys
sys.path.insert(0, "../crisporWebsite")
from crisporEffScores import *

setBinDir("../crisporWebsite/bin")
setCacheDir("./out/")

for fname in glob.glob("effData/*.context.tab"):
    outFname = fname.replace(".context.tab", ".scores.tab")
    if isfile(outFname):
        print "already there, %s" % outFname
        continue
    dataset = basename(fname).split(".")[0]
    print "Processing %s" % fname

    newRows = []
    seqs = []
    for row in iterTsvRows(fname):
        seqs.append(row.longSeq)
        newRow = [dataset, row.guide, row.seq, row.modFreq, row.db, row.pos, row.longSeq]
        newRows.append(newRow)

    scores = calcAllScores(seqs, doAll=True) # removed for the Doench2016 dataset, does not go through the Fusi API
    #scores = calcAllScores(seqs, addOpt=["wangOrig"])
    # add a drsc score to stay compatible with older scripts call housden drsc
    scores["drsc"] = scores["housden"]
    scoreNames = sorted(scores.keys())

    headers = ["dataset", "guide", "seq", "modFreq", "db", "position", "longSeq100Bp"]
    headers.extend(scoreNames)

    ofh = open(outFname, "w")
    writeRow(ofh, headers)
    for i in range(0, len(newRows)):
        row = newRows[i]
        for scoreName in scoreNames:
            scoreList = scores[scoreName]
            row.append(scoreList[i])
        writeRow(ofh, row)
    ofh.close()
    print "wrote %s" % ofh.name
