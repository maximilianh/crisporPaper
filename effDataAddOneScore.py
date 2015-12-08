# add a single score to the effData/*.scores.tab files
# this works by writing all new files to effData.tmp/
# you have then have to move them over manually to effData/

from annotateOffs import *
import glob
import sys
from os.path import *

sys.path.insert(0, "../crispor")
from crisporEffScores import *

# crisprEffScores needs to find its binaries
setBinDir("../crispor/bin")
setCacheDir("./out/")

if not isdir("effData.tmp"):
    os.mkdir("effData.tmp")

tmpFnames = []
fnames = glob.glob("effData/*.scores.tab")
#fnames = ["effData/schoenig.scores.tab"]
for fname in fnames:
    print "reading %s" % fname
    outFname = fname.replace("effData/", "effData.tmp/")
    #if isfile(outFname):
        #print "already there, %s" % outFname
        #continue
    dataset = basename(fname).split(".")[0]
    print "Processing %s" % fname

    newRows = []
    seqs = []
    headers = None
    for row in iterTsvRows(fname):
        headers = list(row._fields)
        seqs.append(row.longSeq100Bp)
        newRow = list(row)
        newRows.append(newRow)

    # CHANGE HERE WHEN ADDING ANOTHER SCORE
    print "Getting scores for %d sequences" % len(seqs)
    shortSeqs = trimSeqs(seqs, -20, 10)
    scoreList = calcWuCrisprScore(shortSeqs)
    #print seqs
    #print scoreList
    #print len(scoreList)
    #print len(seqs)
    headers.append("wuCrispr")
    # END CHANGE

    assert(len(scoreList)==len(newRows))

    ofh = open(outFname, "w")
    writeRow(ofh, headers)
    for i in range(0, len(newRows)):
        row = newRows[i]
        row.append(scoreList[i])
        writeRow(ofh, row)
    ofh.close()
    print "wrote %s" % ofh.name
