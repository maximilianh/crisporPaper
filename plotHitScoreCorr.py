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
    return seqs

def plotScat(ax, otData, onlyStudies=None, onlyRank=False, suppScore=False):
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
            if otSeq in skipSeqs:
                continue

            if suppScore:
                hitScore = calcHsuSuppScore(guideSeq[1:-3], otSeq[1:-3])
            else:
                hitScore = calcHitScore(guideSeq, otSeq)
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
    ax.set_xlabel("Predicted cleavage (rank)")
    ax.set_ylabel("Observed cleavage (rank)")
    #ax.xlim( (0,20) )
    #ax.ylim( (0,0.4) )
    if onlyRank:
        ax.set_xlim((0, 50))
        ax.set_ylim((0, 50))

    # x = scipy.array([-0.65499887,  2.34644428, 3.0])
    # y = scipy.array([-1.46049758,  3.86537321, 21.0])
    r_row, pearP = pearsonr(xVals, yVals)
    rho, spearP = spearmanr(xVals, yVals)
    print "Pearson: r %f, p-Val %f Spearman: rho %f, p-Val %f" % (r_row, pearP, rho, spearP)

def main():

    #guideOffSeqs = parseOfftargetsWithNames("hsu2013/hsuSingle.tab", 10, False, None)
    #plotScat(guideOffSeqs, "out/hitScoreCorr-hsuSingle.pdf")

    fig, axArr = plt.subplots(4,1)
    guideOffSeqs = parseOfftargetsWithNames("out/annotFiltOfftargets.tsv", 10, False, None)
    onlyRank = True
    axArr[0].set_title("Hsu only, MIT score")

    global skipSeqs
    skipSeqs = set(parseHsuSkipped())

    plotScat(axArr[0], guideOffSeqs, onlyStudies=["Hsu"], onlyRank=onlyRank)

    axArr[1].set_title("Hsu only, Hsu Supp score")
    plotScat(axArr[1], guideOffSeqs, onlyStudies=["Hsu"], onlyRank=onlyRank, suppScore=True)

    axArr[2].set_title("All off-targets, MIT score")
    plotScat(axArr[2], guideOffSeqs, onlyStudies=None, onlyRank=onlyRank)

    axArr[3].set_title("All off-targets, Hsu Supp score")
    plotScat(axArr[3], guideOffSeqs, onlyStudies=None, onlyRank=onlyRank, suppScore=True)

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
