# get correlation between the predicted scores for all combinations of scores

from annotateOffs import *
import numpy as np
from itertools import combinations
from scipy.stats import linregress, pearsonr, spearmanr
from os.path import basename

scoreTypes = ["doench", "ssc", "svm", "chariRaw", "finalGc6", "oof", "myScore"]

def parseScores(inDir, datasets):
    " index all scores in inDir and return as dataset -> scoreType -> list of scores "
    # first get all seqs and calc all seq scores
    print "parsing file in %s" % inDir
    seqKoEffs = {}
    for dataset in datasets:
        fname = "effData/%s.ext.tab" % dataset
        #for fname in glob.glob(inDir+"/*.ext.tab"):
            #dataset = basename(fname).split(".")[0]
        #if dataset not in datasets:
            #continue
        if "doench2014-Mm" in fname:
            continue # lots of duplicates in this dataset
        for row in iterTsvRows(fname):
            freq = float(row.modFreq)
            seqKoEffs[row.extSeq.upper()] = freq

    seqScores = calcEffScores(seqKoEffs)

    # now split them into scoreType -> list of scores
    guideData = defaultdict(list)
    scoreTypes = seqScores.values()[0]
    for seq, scores in seqScores.iteritems():
        for scoreType in scoreTypes:
            score = scores[scoreType]
            guideData[scoreType].append(score)
    return guideData

    # now split them into dataset -> scoreType -> list of scores
    #fileGuideData = {}
    #fileMinReqFreqs = {}
    #for dataset in datasets:
        #fname = "effData/%s.ext.tab" % dataset
        ##for fname in glob.glob(inDir+"/*.ext.tab"):
        ##dataset = basename(fname).split(".")[0]
        ##if dataset not in datasets:
            ##continue
        #guideData = defaultdict(list)
        #koEffs = []
        #for row in iterTsvRows(fname):
            #seq = row.extSeq
            #scores = seqScores[seq]
            #for scoreType, score in scores.iteritems():
                #guideData[scoreType].append( score )
#
        #dataName = basename(fname).split(".")[0]
        #fileGuideData[dataName] = guideData
    #return fileGuideData

def main():
    scoreDict = parseScores("effData", ["xu2015Train", "chari2015Train", "doench2014-Hs"])

    rows = []
    #for dataset, scoreDict in dataScores.iteritems():
        #if "xu2015Train" not in dataset:
            #continue
    for score1, score2 in combinations(scoreTypes, 2):
        scores1 = scoreDict[score1]
        scores2 = scoreDict[score2]
        if None in scores1 or None in scores2:
            print "misssing scores, skipping %s-%s" % (score1, score2)
            continue
        #print scores1
        #print scores2
        pearR, pVal = pearsonr(scores1, scores2)
        spearR, pVal = spearmanr(scores1, scores2)

        #print dataset, score1, score2, pearR, spearR
        row = [score1, score2, "%0.3f" % pearR, "%0.3f" % spearR]
        rows.append(row)
        #print row

    rows.sort(key=operator.itemgetter(3), reverse=True)

    ofh = open('out/scoreCorr.tsv', "w")
    ofh.write("score1\tscore2\tpearsonCorr\tspearmanCorr\n")
    for r in rows:
        ofh.write( "\t".join(r))
        ofh.write("\n")
    print "written to %s" % ofh.name

print "Chari et al report a correlation of 0.819-853 Fig2a"
main()
