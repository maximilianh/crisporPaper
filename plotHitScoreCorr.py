# plot the correlation between the MIT hitScore and the cleavage frequency on the Hsu et al 2013
# data or our meta dataset

import os, logging, itertools
logging.basicConfig(loglevel=logging.INFO)
from os.path import isfile, splitext, join
from annotateOffs import *
from collections import defaultdict

logging.basicConfig(loglevel=logging.INFO)

from scipy.stats import linregress
from scipy.stats import pearsonr, spearmanr
import matplotlib.pyplot as plt

import crisporOtScores

# add guide names to scatter dots?
doLabels = False

skipSeqs = None

def parseHsuSkipped():
    " return a list of the Hsu off-targets that have < 10000 reads "
    seqs= []
    for line in open("origData/Hsu2013/readCounts.tab"):
        if line.startswith("otSeq"):
            continue
        seq, readCount = line.strip().split()
        readCount = int(readCount)
        if readCount < 10000:
            seqs.append(seq)
    print "Found %d Hsu seqs that will be skipped" % len(seqs)
    return seqs

def plotScat(ax, otData, onlyStudies=None, onlyRank=False, scoreName="mit"):
    xVals, yVals = [], []
    figs = []
    markers = itertools.cycle(["*", "+", "o", "x", "^"])
    colors = itertools.cycle(["black", "green", "red", "blue"])
    guideNames = []

    i = 0
    for (guideName, guideSeq), otSeqs in otData.iteritems():
        study = guideName.split("_")[0].split("/")[0]
        if onlyStudies is not None and study not in onlyStudies:
            #print "skipping guide", guideName
            continue
        #print "processing guide", guideName
        for otSeq, otFreq in otSeqs.iteritems():
            if "-" in otSeq:
                #print "OT-Gap: skipping %s / %s " % (guideName, otSeq)
                continue

            if otSeq in skipSeqs:
                continue

            if scoreName=="hsuSupp":
                #hitScore = calcHsuSuppScore(guideSeq[1:-3], otSeq[1:-3])
                hitScore = crisporOtScores.calcHsuSuppScore(guideSeq, otSeq, strat="col")
                #print "running Hsu with %s=%s and %s" % (guideName, guideSeq, otSeq)
                #hitScore = crisporOtScores.calcRawHsu(guideSeq, otSeq)
                #hitScore = crisporOtScores.calcHsuSuppScore2(guideSeq, otSeq)
            elif scoreName=="mit":
                hitScore = calcHitScore(guideSeq, otSeq)
            elif scoreName=="cfd":
                hitScore = calcCfdScore(guideSeq, otSeq)
            else:
                print scoreName
                assert(False)
            #print hitScore, guideSeq, otSeq
            xVals.append(hitScore)
            yVals.append(otFreq)

        if onlyRank:
            xVals = useRanks(xVals)
            yVals = useRanks(yVals)

        if doLabels:
            for hitScore, otFreq in zip(xVals, yVals):
                annots.append( (hitScore, otFreq, guideName) )

        marker = markers.next()
        color = colors.next()
        fig = ax.scatter(xVals, yVals, \
            alpha=.5, \
            marker = marker, \
            color = color, \
            s=5)
        figs.append(fig)
        #doLabels.append(title)
        i+=1

    if doLabels:
        for x, y, label in annots:
           ax.annotate(
              label, fontsize=7, rotation=0, ha="right", rotation_mode="anchor",
              xy = (x, y), xytext = (0,0), alpha=0.9,
              textcoords = 'offset points', va = 'bottom')
    #ax.legend(figs,
           #studyNames,
           #scatterpoints=1,
           #loc='upper right',
           #ncol=3,
           #fontsize=10)
    if onlyRank:
        ax.set_xlabel("Predicted cleavage (rank)")
        ax.set_ylabel("Observed cleavage (rank)")
        ax.set_xlim((0, 50))
        ax.set_ylim((0, 50))
    else:
        ax.set_xlabel("Predicted cleavage score")
        ax.set_ylabel("Modif. frequency")
        #ax.set_xlim( (0,1.0) )
        #ax.set_ylim( (0,0.3) )

    # x = scipy.array([-0.65499887,  2.34644428, 3.0])
    # y = scipy.array([-1.46049758,  3.86537321, 21.0])
    r_row, pearP = pearsonr(xVals, yVals)
    rho, spearP = spearmanr(xVals, yVals)
    print "score %s: Pearson: r %f, p-Val %f Spearman: rho %f, p-Val %f" % (scoreName, r_row, pearP, rho, spearP)

def main():

    #guideOffSeqs = parseOfftargetsWithNames("hsu2013/hsuSingle.tab", 10, False, None)
    #plotScat(guideOffSeqs, "out/hitScoreCorr-hsuSingle.pdf")

    fig, axArr = plt.subplots(5,1)
    fname = "out/annotOfftargets.tsv"
    #fname = "out/annotOfftargets.tsv"
    guideOffSeqs = parseOfftargetsWithNames(fname, 10, False, None)
    onlyRank = True
    axArr[0].set_title("Hsu only, MIT score")

    global skipSeqs
    skipSeqs = set(parseHsuSkipped())

    #plotScat(axArr[0], guideOffSeqs, onlyStudies=["Hsu"], onlyRank=onlyRank)

    axArr[0].set_title("Hsu only, Hsu Supp score")
    plotScat(axArr[0], guideOffSeqs, onlyStudies=["Hsu"], onlyRank=onlyRank, scoreName="hsuSupp")

    axArr[1].set_title("Hsu only, MIT score")
    plotScat(axArr[1], guideOffSeqs, onlyStudies=["Hsu"], onlyRank=onlyRank, scoreName="mit")

    axArr[2].set_title("Hsu only, Cfd score")
    plotScat(axArr[2], guideOffSeqs, onlyStudies=["Hsu"], onlyRank=onlyRank, scoreName="cfd")

    #axArr[2].set_title("All off-targets, Hsu Supp score")
    #plotScat(axArr[2], guideOffSeqs, onlyStudies=None, onlyRank=onlyRank, scoreName="hsuSupp")

    #axArr[3].set_title("All off-targets, MIT score")
    #plotScat(axArr[3], guideOffSeqs, onlyStudies=None, onlyRank=onlyRank, scoreName="mit")

    #axArr[4].set_title("All off-targets, CFD score")
    #plotScat(axArr[4], guideOffSeqs, onlyStudies=None, onlyRank=onlyRank, scoreName="cfd")

    outPdfName = "out/hitScoreCorr-allFilt.pdf"
    fig.set_size_inches(5,15)
    fig.tight_layout()
    plt.savefig(outPdfName)
    plt.savefig(outPdfName.replace("pdf", "png"))
    print "plot written to %s and .png" % outPdfName
    plt.close()

    #for (guideName, guideSeq), otSeqs in guideOffSeqs.iteritems():
        #for otSeq, otFreq in otSeqs.iteritems():
            #hitScore = calcHitScore(guideSeq, otSeq)
            #print guideName, guideSeq, hitScore, otFreq

main()
