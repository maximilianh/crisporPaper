# get correlation between the predicted scores for all combinations of scores

from annotateOffs import *
import numpy as np
from itertools import combinations
from scipy.stats import linregress, pearsonr, spearmanr

scoreTypes = ["doench", "ssc", "svm", "chariRaw", "finalGc6"]

def parseScores(inDir):
    " index all scores in inDir and return as dataset -> scoreType -> list of scores "
    # first get all seqs and calc all seq scores
    print "parsing file in %s" % inDir
    seqKoEffs = {}
    for fname in glob.glob(inDir+"/*.ext.tab"):
        for row in iterTsvRows(fname):
            freq = float(row.modFreq)
            seqKoEffs[row.extSeq] = freq

    seqScores = calcEffScores(seqKoEffs)

    # now split them into dataset -> scoreType -> list of scores
    fileGuideData = {}
    fileMinReqFreqs = {}
    for fname in glob.glob(inDir+"/*.ext.tab"):
        guideData = defaultdict(list)
        koEffs = []
        for row in iterTsvRows(fname):
            seq = row.extSeq
            scores = seqScores[seq]
            for scoreType, score in scores.iteritems():
                guideData[scoreType].append( score )

        dataName = basename(fname).split(".")[0]
        fileGuideData[dataName] = guideData
    return fileGuideData

def main():
    dataScores = parseScores("effData")

    rows = []
    for dataset, scoreDict in dataScores.iteritems():
        if "xu2015Train" not in dataset:
            continue
        for score1, score2 in combinations(scoreTypes, 2):
            scores1 = scoreDict[score1]
            scores2 = scoreDict[score2]
            pearR, pVal = pearsonr(scores1, scores2)
            spearR, pVal = spearmanr(scores1, scores2)

            row = [dataset, score1, score2, "%0.3f" % pearR, "%0.3f" % spearR]
            rows.append(row)
            #print row

    rows.sort(key=operator.itemgetter(3), reverse=True)

    ofh = open('out/scoreCorr.tsv', "w")
    ofh.write("dataset\tscore1\tscore2\tpearsonCorr\tspearmanCorr\n")
    for r in rows:
        ofh.write( "\t".join(r))
        ofh.write("\n")
    print "written to %s" % ofh.name

main()
