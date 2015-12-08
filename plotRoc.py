# plot the ROC curve of the off-targets when varying the off-target score
import glob, sys, logging
from annotateOffs import *

import numpy as np
import matplotlib.pyplot as plt
from os.path import basename, splitext, isfile
from sklearn.metrics import roc_curve, auc

import crisporOtScores

# we only look at off-targets with a certain number of mismatches
# otherwise the ROC curve would not go up to 1.0 as MIT can only search for 4 MMs
maxMismatches = 4

# for the ROC curve, we only analyze off-targets with certain PAM sites 
# assuming that no software can find the relatively rare PAM sites
# that are not GG/GA/AG
validPams = ["GG", "GA", "AG"]
#validPams = ["GG"]

# only look at alternative PAMs, can be used to determine best cutoff for the alternative PAMs
#onlyAlt = True
onlyAlt = False

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

def plotRoc(prefix, allOts, expFreqs, predScores, colors, styles, plots, labels, isCropit=False, fracs=[0.0], legLabels=[""]):
    " plot ROC curve and write annotation to ofh file "
    print "plotting ROC for %s" % prefix
    i= 0
    maxSens = 0
    #fracList = [0.01, 0.001]
    #fracList = [0.00]

    for label, minFrac in zip(legLabels, fracs):
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
        print "AUC=%f" % roc_auc
        if label!="":
            plotLabel = label + ", AUC %0.2f" % roc_auc
        elif minFrac == 0.0:
            plotLabel = prefix + ", AUC %0.2f" % roc_auc
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

def compHsuMit(guideValidOts, guideSeqs, outFname):
    " calc MIT and hsu OT scores and write to outfname "
    ofh = open(outFname, "w")
    ofh.write("\t".join(["guide", "offt", "modFreq", "MIT", "hsu"])+"\n")
    for guideSeq, otFreqs in guideValidOts.iteritems():
        for otSeq, otFreq in otFreqs.iteritems():
            mitScore = crisporOtScores.calcMitScore(guideSeq, otSeq)
            hsuScore = crisporOtScores.calcHsuSuppScore2(guideSeq, otSeq)
            row = [guideSeq, otSeq, otFreq, mitScore, hsuScore]
            row = [str(x) for x in row]
            ofh.write("\t".join(row)+"\n")
    print "wrote %s" % ofh.name

        
def main():
    guideValidOts, guideSeqs = parseOfftargets("out/annotFiltOfftargets.tsv", maxMismatches, onlyAlt, validPams)

    tmpFname = "/tmp/crisporOffs.pickle"

    mitPredOts = parseMit("mitOfftargets", guideSeqs)

    compHsuMit(guideValidOts, guideSeqs, "out/hsuVsMit.tsv")

    if isfile(tmpFname):
        print "reading offtargets from %s" % tmpFname
        crisporPredOts = pickle.load(open(tmpFname))
    else:
        crisporPredOts = parseCrispor("crisporOfftargets", guideSeqs, maxMismatches)
        pickle.dump(crisporPredOts, open(tmpFname, "w"))
        print "Wrote offtargets to %s" % tmpFname

    #cropitPredOts = parseCropit("cropitOfftargets", guideSeqs)

    plots = []
    labels = []

    plt.figure(figsize=(7,7))
    #dataName = "filtered BWA"

    colors = ["darkblue", "red"]
    styles = ["-", "-", ":"]
    dataName = "MIT score"
    crisporScores = collapseDicts(crisporPredOts)
    validOffts = collapseDicts(guideValidOts)
    plots, labels = plotRoc(dataName, crisporScores, validOffts, crisporScores, colors, styles, plots, labels, fracs=[0.0, 0.01], legLabels=["MIT score", "MIT score (freq. > 1%)"])
    if not onlyAlt:
        #plt.annotate('True Pos. Rate = 0.96:\nMIT score = 0.1', xy=(0.588221664414, 0.9640), xytext=(0.4, .84),
        plt.annotate('True Pos. Rate = 1.0 / 0.96:\nMIT score = 0.1', xy=(0.588221664414, 0.9650), xytext=(0.5, 1.04),
                arrowprops=dict(facecolor='black', arrowstyle="->"), annotation_clip=False)
        plt.annotate('True Pos. Rate = 1.0 / 0.96:\nMIT score = 0.1', xy=(0.588221664414, 1.0), xytext=(0.5, 1.04),
                arrowprops=dict(facecolor='black', arrowstyle="->"), annotation_clip=False)

    colors = ["blue", "green"]
    styles = ["-", ":", ":"]
    mitScores = collapseDicts(mitPredOts)
    plots, labels = plotRoc("MIT Website", crisporScores, validOffts, mitScores, colors, styles, plots, labels)

    #plt.annotate('offt. score = ?', xy=(0.8, 1.0), xytext=(0.7, .84),

                #arrowprops=dict(facecolor='black', arrowstyle="->"))

    #colors = ["black"]
    #styles = ["-", "-", "-"]
    #cfdPredOts = calcOtScores(crisporPredOts, calcCfdScore)
    #cfdScores = collapseDicts(cfdPredOts)
    #plots, labels = plotRoc("CFD score", cfdScores, validOffts, cfdScores, colors, styles, plots, labels)

    cropitScores = collapseDicts(calcOtScores(crisporPredOts, calcCropitScore))
    styles = ["-.", "-."]
    colors = ["red"]
    plots, labels = plotRoc("Cropit score", cropitScores, validOffts, cropitScores, colors, styles, plots, labels)

    ccTopPredOts = calcOtScores(crisporPredOts, calcCcTopScore)
    ccTopScores = collapseDicts(ccTopPredOts)
    styles = ["--", "--"]
    colors = ["green"]
    plots, labels = plotRoc("CCTop score", ccTopScores, validOffts, ccTopScores, colors, styles, plots, labels)


    hsuScores = collapseDicts(calcOtScores(crisporPredOts, crisporOtScores.calcHsuSuppScore2))
    styles = ["--", "--"]
    colors = ["grey"]
    plots, labels = plotRoc("Hsu score", hsuScores, validOffts, hsuScores, colors, styles, plots, labels)

    #for strat in ["raw", "all", "none", "row", "col", "onlyAvgs", "avgs", "none_Sum", "none_allSum", "limit"]:
        #print "Hsu normalisation: %s" % strat
        #crisporOtScores.loadHsuMat(strat)
        #hsuScores = collapseDicts(calcOtScores(crisporPredOts, crisporOtScores.calcHsuSuppScore))
        #styles = ["-", "-"]
        #plots, labels = plotRoc("Hsu score", hsuScores, validOffts, hsuScores, colors, styles, plots, labels)


    #colors = ["black", "blue", "green"]
    #styles = ["--", "--", "--"]
    #plots, labels = plotRoc("CROP-IT", crisporScores, validOffts, cropitPredOts, colors, styles, plots, labels, isCropit=True)

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
