# plot the correlation between the MIT hitScore and the cleavage frequency on both the Hsu et al 2013
# data and our meta dataset

import os, logging, itertools
logging.basicConfig(loglevel=logging.INFO)
from os.path import isfile, splitext, join
from annotateOffs import *
from collections import defaultdict

logging.basicConfig(loglevel=logging.INFO)

from scipy.stats import linregress
from scipy.stats import pearsonr
import matplotlib.pyplot as plt

def plotScat(otData, outPdfName):
    xVals, yVals = [], []
    figs = []
    annots = []
    markers = itertools.cycle(["*", "+", "o", "x", "-", "^"])
    colors = itertools.cycle(["black", "green", "red", "blue"])

    for (guideName, guideSeq), otSeqs in otData.iteritems():
        for otSeq, otFreq in otSeqs.iteritems():
            hitScore = calcHitScore(guideSeq, otSeq)
            xVals.append(hitScore)
            yVals.append(otFreq)
            #annots.append( (hitScore, otFreq, guideName) )

        fig = plt.scatter(xVals, yVals, \
            alpha=.5, \
            marker = markers[i], \
            color = colors[i], \
            s=5)
        figs.append(fig)
        #labels.append(title)

    for x, y, label in annots:
       plt.annotate(
          label, fontsize=7, rotation=30, ha="right", rotation_mode="anchor",
          xy = (x, y), xytext = (0,0), alpha=0.9,
          textcoords = 'offset points', va = 'bottom')
    plt.xlabel("Hit Score")
    plt.ylabel("Modification Frequency")
    #plt.xlim( (0,20) )
    #plt.ylim( (0,0.4) )
    plt.savefig(outPdfName)
    plt.savefig(outPdfName.replace("pdf", "png"))
    print "plot written to %s and .png" % outPdfName

    # x = scipy.array([-0.65499887,  2.34644428, 3.0])
    # y = scipy.array([-1.46049758,  3.86537321, 21.0])
    r_row, p_value = pearsonr(xVals, yVals)
    print "r, p", r_row, p_value

def main():
    #guideOffSeqs = parseOfftargetsWithNames("out/annotFiltOfftargets.tsv", 10, False, ["GG","AG","GA"])
    guideOffSeqs = parseOfftargetsWithNames("hsu2013/hsuSingle.tab", 10, False, None)
    plotScat(guideOffSeqs, "out/hitScoreCorr-hsuSingle.pdf")

    guideOffSeqs = parseOfftargetsWithNames("out/annotFiltOfftargets.tsv", 10, False, None)
    plotScat(guideOffSeqs, "out/hitScoreCorr-allFilt.pdf")

    #for (guideName, guideSeq), otSeqs in guideOffSeqs.iteritems():
        #for otSeq, otFreq in otSeqs.iteritems():
            #hitScore = calcHitScore(guideSeq, otSeq)
            #print guideName, guideSeq, hitScore, otFreq

main()
