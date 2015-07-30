import os, logging
from annotateOffs import *
from collections import defaultdict
#from scipy.stats import linregress

import matplotlib.pyplot as plt
import numpy as np

def parseAllSeqs(dirName):
    """ return a list of all 34-mer sequences for which we have scores
    and also a list of (dataname, list of KO efficiencies)
    """
    seqs = set()
    effLists = []
    for fname in glob.glob(dirName+"/*.ext.tab"):
        if "S7" in fname or "S10" in fname:
            continue
        effList =[]
        for row in iterTsvRows(fname):
            seqs.add(row.extSeq)
            effList.append(float(row.modFreq))
        effLists.append(( basename(fname).split(".")[0] , effList) )
    return seqs, effLists

def main():
    #plt.figure(figsize=(,10))
    #fig, axArr = plt.subplots(4, 1, sharex="col")
    fig, (axRow1, axRow2) = plt.subplots(2, 6)
    fig.set_size_inches(20,5)

    seqs, effLists = parseAllSeqs("effData")
    scores = calcEffScores(seqs)
    
    scoresByType = defaultdict(list)
    for seq, seqScores in scores.iteritems():
        for seqType, score in seqScores.iteritems():
            scoresByType[seqType].append(score)

    # scores
    for plotRow, scoreType in enumerate(["svm", "doench", "ssc", "chariRaw", "finalGc6"]):
            seqScores = scoresByType[scoreType]

            ax = axRow1[plotRow]
            ax.hist(seqScores)
            ax.set_xlabel("%s score" % scoreType)
            ax.set_ylabel("Frequency")
            #ax.set_ylim(0,1.0)
            #if plotCol==0:
                #ax.set_title(gene)

    # efficiencies
    for plotRow, (dataName, effList) in enumerate(effLists):
        ax = axRow2[plotRow]
        ax.hist(effList)
        ax.set_xlabel("KO efficiency")
        ax.set_ylabel("Frequency")
        ax.set_title(dataName)

    fig.tight_layout()

    outFname = "out/effScoreHistogram.pdf"
    plt.savefig(outFname)
    plt.savefig(outFname.replace(".pdf", ".png"))
    print "wrote plot to %s and .png" % outFname

main()
