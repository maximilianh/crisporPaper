# plot the ROC curve of the off-targets when varying the off-target score
import glob, sys, logging
from annotateOffs import *

import numpy as np
import matplotlib.pyplot as plt
from os.path import basename, splitext, isfile
from sklearn.metrics import roc_curve, auc

# we only look at off-targets with a certain number of mismatches
# otherwise the ROC curve would not go up to 1.0 as MIT can only search for 4 MMs
maxMismatches = 4

# for the ROC curve, we only analyze off-targets with certain PAM sites 
# assuming that no software can find the relatively rare PAM sites
# that are not GG/GA/AG
validPams = ["GG", "GA", "AG"]

# only look at alternative PAMs, can be used to determine best cutoff for the alternative PAMs
onlyAlt = True

if len(sys.argv)>1:
    altPamCutoff = float(sys.argv[1])
else:
    altPamCutoff = None

def parseCropit(inDir, guideSeqs):
    " parse the cropit minimal files, return a dict with guideSeq -> otSeq -> otScore "
    data = defaultdict(dict)
    for guideName in guideSeqs:
        guideNameNoCell = guideName.replace("/K562", "").replace("/Hap1","")
        fname = join(inDir, guideNameNoCell+".tsv")
        print "parsing %s" % fname
        if not isfile(fname):
            logging.error("MISSING: %s" % fname)
            continue
        for line in open(fname):
            fs = line.strip("\n").split('\t')
            otSeq = fs[0]
            score = fs[1]
            guideSeq = guideSeqs[guideName]
            data[guideSeq][otSeq]=float(score)
    return data

def plotRoc(prefix, allOts, expFreqs, predScores, colors, styles, plots, labels, isCropit=False):
    " plot ROC curve and write annotation to ofh file "
    print "plotting ROC for %s" % prefix
    i= 0
    maxSens = 0
    fracList = [0.01, 0.001]

    for minFrac in fracList:
        yVals = []
        yTrue = []
        #for guideSeq, guideSeqFreqs in guideValidOts.iteritems():
        for otSeq in allOts:
            otScore = predScores.get(otSeq, 0.0)
            yVals.append(otScore)
            if expFreqs.get(otSeq, 0.0) > minFrac:
                yTrue.append(1)
            else:
                yTrue.append(0)

        #print yVals
        #print yTrue
        assert(len(yTrue)==len(yVals))
        print "%s, minFreq %f, set of elements: %d "% (prefix, minFrac, len(yVals))
        print("number of validated off-targets: %d" % (sum(yTrue)))
        fpr, tpr, thresholds = roc_curve(yTrue, yVals)
        if prefix == "CRISPOR":
            for f, t, score in zip(fpr, tpr, thresholds):
                print minFrac, f, t, score
        roc_auc = auc(fpr, tpr)
        if minFrac == 0.0:
            plotLabel = prefix+", no freq. limit (%d off-targets)" % (len(yTrue))
        else:
            plotLabel = prefix+", mod. freq. > %0.1f%%: AUC %0.2f" % ((minFrac*100), roc_auc)
        p, = plt.plot(fpr, tpr, ls=styles[i], color=colors[i]) # NB: keep the comma!
        plots.append(p)
        labels.append(plotLabel)
        i+=1
    return plots, labels

def collapseDicts(predOts):
    " given a nested dict guideSeq -> otSeq -> score, return a merged otSeq -> score "
    ret = {}
    for guideSeq, predOtScores in predOts.iteritems():
        for otSeq, otScore in predOtScores.iteritems():
            ret[otSeq] = otScore
    return ret

def main():
    guideValidOts, guideSeqs = parseOfftargets("out/annotFiltOfftargets.tsv", maxMismatches, onlyAlt, validPams)
    crisporPredOts = parseCrispor("crisporOfftargets", guideSeqs, maxMismatches)
    mitPredOts = parseMit("mitOfftargets", guideSeqs)
    cropitPredOts = parseCropit("cropitOfftargets", guideSeqs)

    plots = []
    labels = []

    colors = ["black", "blue", "green"]
    styles = ["-", "-", "-"]

    plt.figure(figsize=(7,7))
    #dataName = "filtered BWA"
    dataName = "CRISPOR"
    crisporScores = collapseDicts(crisporPredOts)
    validOffts = collapseDicts(guideValidOts)
    plots, labels = plotRoc(dataName, crisporScores, validOffts, crisporScores, colors, styles, plots, labels)
    if not onlyAlt:
        plt.annotate('TPR = 0.96:\nofft. score = 0.1', xy=(0.588221664414, 0.9640), xytext=(0.4, .84),
                arrowprops=dict(facecolor='black', arrowstyle="->"))
    #plt.annotate('offt. score = ?', xy=(0.8, 1.0), xytext=(0.7, .84),
                #arrowprops=dict(facecolor='black', arrowstyle="->"))

    colors = ["black", "blue", "green"]
    styles = [":", ":", ":"]
    mitScores = collapseDicts(mitPredOts)
    plots, labels = plotRoc("MIT", crisporScores, validOffts, mitScores, colors, styles, plots, labels)

    #colors = ["black", "blue", "green"]
    #styles = ["--", "--", "--"]
    plots, labels = plotRoc("CROP-IT", crisporScores, validOffts, cropitPredOts, colors, styles, plots, labels, isCropit=True)

    plt.legend(plots,
           labels,
           loc='lower right',
           ncol=1,
           fontsize=11)

    plt.xlabel("False positive rate")
    plt.ylabel("True positive rate")

    #ax = plt.gca()
    #ax.axhline(y=maxSens1, ls=":", color="k")
    #ax.axvline(x=0.6, ls="-", color="b")

    #plt.text(0, maxSens1, "max = %0.2f" % maxSens1)
    plt.ylim(0,1.0)
    plt.xlim(0,1.0)

    outfname = "out/roc.pdf"
    plt.savefig(outfname)
    print "wrote %s" % outfname

    outfname = "out/roc.png"
    plt.savefig(outfname)
    print "wrote %s" % outfname
    #print "wrote data to %s" % ofh.name

main()
