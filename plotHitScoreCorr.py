# plot the correlation between the MIT hitScore and the cleavage frequency on both the Hsu et al 2013
# data and our meta dataset

import os, logging, itertools
logging.basicConfig(loglevel=logging.INFO)
from os.path import isfile, splitext, join
from annotateOffs import *
from collections import defaultdict

logging.basicConfig(loglevel=logging.INFO)

from scipy.stats import linregress
from scipy.stats import pearsonr, spearmanr
import matplotlib.pyplot as plt

def plotScat(otData, outPdfName, onlyStudies=None, onlyRank=False):
    xVals, yVals = [], []
    figs = []
    annots = []
    markers = itertools.cycle(["*", "+", "o", "x", "^"])
    colors = itertools.cycle(["black", "green", "red", "blue"])
    guideNames = []

    i = 0
    for (guideName, guideSeq), otSeqs in otData.iteritems():
        study = guideName.split("_")[0].split("/")[0]
        if onlyStudies is not None and study not in onlyStudies:
            #print "skipping guide", guideName
            continue
        print "processing guide", guideName
        for otSeq, otFreq in otSeqs.iteritems():
            hitScore = calcHitScore(guideSeq, otSeq)
            xVals.append(hitScore)
            yVals.append(otFreq)

        if onlyRank:
            xVals = useRanks(xVals)
            yVals = useRanks(yVals)

        #for hitScore, otFreq in zip(xVals, yVals):
            #annots.append( (hitScore, otFreq, guideName) )
        marker = markers.next()
        color = colors.next()
        fig = plt.scatter(xVals, yVals, \
            alpha=.5, \
            marker = marker, \
            color = color, \
            s=5)
        figs.append(fig)
        #labels.append(title)
        i+=1

    for x, y, label in annots:
       plt.annotate(
          label, fontsize=7, rotation=30, ha="right", rotation_mode="anchor",
          xy = (x, y), xytext = (0,0), alpha=0.9,
          textcoords = 'offset points', va = 'bottom')
    #plt.legend(figs,
           #studyNames,
           #scatterpoints=1,
           #loc='upper right',
           #ncol=3,
           #fontsize=10)
    plt.xlabel("Predicted cleavage (rank)")
    plt.ylabel("Observed cleavage (rank)")
    #plt.xlim( (0,20) )
    #plt.ylim( (0,0.4) )
    plt.savefig(outPdfName)
    plt.savefig(outPdfName.replace("pdf", "png"))
    print "plot written to %s and .png" % outPdfName
    plt.close()

    # x = scipy.array([-0.65499887,  2.34644428, 3.0])
    # y = scipy.array([-1.46049758,  3.86537321, 21.0])
    r_row, pearP = pearsonr(xVals, yVals)
    rho, spearP = spearmanr(xVals, yVals)
    print "Pearson: r %f, p-Val %f Spearman: rho %f, p-Val %f" % (r_row, pearP, rho, spearP)

def main():

    #guideOffSeqs = parseOfftargetsWithNames("hsu2013/hsuSingle.tab", 10, False, None)
    #plotScat(guideOffSeqs, "out/hitScoreCorr-hsuSingle.pdf")

    guideOffSeqs = parseOfftargetsWithNames("out/annotFiltOfftargets.tsv", 10, False, None)
    plotScat(guideOffSeqs, "out/hitScoreCorr-allFilt.pdf", onlyStudies=["Hsu"], onlyRank=True)

    #for (guideName, guideSeq), otSeqs in guideOffSeqs.iteritems():
        #for otSeq, otFreq in otSeqs.iteritems():
            #hitScore = calcHitScore(guideSeq, otSeq)
            #print guideName, guideSeq, hitScore, otFreq

main()
