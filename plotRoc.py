import glob, sys, logging
from annotateOffs import *

import numpy as np
import matplotlib.pyplot as plt
from os.path import basename, splitext, isfile

# we only look at off-targets with a certain number of mismatches
maxMismatches = 4

# for the ROC curve, we only analyze off-targets with certain PAM sites 
# assuming that no software can find the relatively rare PAM sites
# that are not GG/GA/AG
validPams = ["GG", "GA", "AG"]

# !!!
# only look at alternative PAMs, can be used to determine best cutoff for the alternative PAMs
# supplemental data??
onlyAlt = False

# the cutoff is varied over these values for the MIT score
mitScoreList = [0.0, 0.001, 0.002, 0.003, 0.005, 0.007, 0.009, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06, 0.07, 0.08, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, 3, 4, 5, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30, 35, 40, 50, 60]

# for the cropit score, we use these values instead
cropitScoreList = range(150, 700, 20)

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

def filterValidOfftargets(guideValidOts, minReadFrac):
    " return a list of all validated offtargets with minReadFrac "
    validOts = set()
    for guideSeq, validOtSeqs in guideValidOts.iteritems():
        for seq, readFrac in validOtSeqs.iteritems():
            if readFrac > minReadFrac:
                validOts.add(seq)
    return validOts

def getRocValues(toolName, guideValidOts, guidePredOts, minReadFrac, ofh, isCropit=False):
    " return a list of (sens, fdr) tuples for a ROC curve plot and output rows to ofh"
    # keep only the validated off-targets with read fraction > minCutoff
    validOts = filterValidOfftargets(guideValidOts, minReadFrac)

    sensList = []
    fdrList = []

    cutoffs = mitScoreList
    if isCropit:
        cutoffs = cropitScoreList

    for cutoff in cutoffs:
        print "XX cutoff", cutoff
        predOts = set()
        allOts = set()
        allOts.update(validOts) # make a copy of the elements

        for guideSeq, predSeqScores in guidePredOts.iteritems():
            # get the predicted sequences over the off-target score cutoff
            for predSeq, seqScore in predSeqScores.iteritems():
                allOts.add(predSeq)

                # check if alternative PAM
                if predSeq[-2:] in ["AG", "GA"] and altPamCutoff!=None and seqScore < altPamCutoff:
                    continue
                elif onlyAlt:
                    continue

                if seqScore > float(cutoff):
                    predOts.add(predSeq)

        notPredOts = allOts - predOts
        if cutoff==0.0 and toolName=="CRISPOR":
            print "missed off-targets by crispor for mod freq > %f: %s" % (minReadFrac, notPredOts)
        notValidOts = allOts - validOts

        tp = validOts.intersection(predOts)
        tn = notPredOts.intersection(notValidOts)
        fp = predOts - validOts
        fn = notPredOts.intersection(validOts)

        # sensitivity - proportion of validated seqs that predicted to be off-targets
        # relative to all off-targets
        sens = float(len(tp)) / (len(tp)+len(fn))

        # specificity - proportion of that are predicted to be not off-targets
        if len(tn)+len(fp)!=0:
            spec = float(len(tn)) / (len(tn)+len(fp))
        else:
            spec = 0.0
        fdr = 1.0 - spec

        sensList.append(sens)
        fdrList.append(fdr)

        row = [toolName, minReadFrac, cutoff, sens*100, fdr, len(tp), len(fp), len(fn), len(tn)]
        row = [str(x) for x in row]
        ofh.write("\t".join(row))
        ofh.write("\n")
        #sys.stdout.flush()

    return sensList, fdrList, validOts

def plotRoc(prefix, guideValidOts, guidePredOts, colors, styles, plots, labels, ofh, isCropit=False):
    " plot ROC curve and write annotation to ofh file "
    i= 0
    maxSens = 0
    fracList = [0.0, 0.001, 0.01]
    if isCropit:
        fracList = [0.01]

    for minFrac in fracList:
        sensList, fdrList, validSeqs = getRocValues(prefix, guideValidOts, guidePredOts, minFrac, ofh, isCropit)
        if minFrac == 0.0:
            plotLabel = prefix+", no freq. limit (%d off-targets)" % (len(validSeqs))
        else:
            plotLabel = prefix+", mod. freq. > %0.1f%% (%d off-targets)" % ((minFrac*100), len(validSeqs))
        p, = plt.plot(fdrList, sensList, ls=styles[i], color=colors[i]) # NB: comma!
        plots.append(p)
        labels.append(plotLabel)
        maxSens = max(maxSens, max(sensList))
        i+=1
    return plots, labels, maxSens

def main():
    guideValidOts, guideSeqs = parseOfftargets("out/annotFiltOfftargets.tsv", maxMismatches, onlyAlt, validPams)
    guidePredOts = parseCrispor("crisporOfftargets", guideSeqs, maxMismatches)
    mitPredOts = parseMit("mitOfftargets", guideSeqs)
    cropitPredOts = parseCropit("cropitOfftargets", guideSeqs)

    ofh = open("out/rocData.tsv", "w")
    headers = ["dataset", "readFrac", "cutoff", "sensitivity", "falseDescRate", "TP", "FP", "FN", "TN"]
    ofh.write("\t".join(headers)+"\n")

    plots = []
    labels = []

    colors = ["black", "blue", "green"]
    styles = ["-", "-", "-"]

    plt.figure(figsize=(7,7))
    #dataName = "filtered BWA"
    dataName = "CRISPOR"
    plots, labels, maxSens1 = plotRoc(dataName, guideValidOts, guidePredOts, colors, styles, plots, labels, ofh)

    colors = ["black", "blue", "green"]
    styles = [":", ":", ":"]
    plots, labels, maxSens2 = plotRoc("MIT", guideValidOts, mitPredOts, colors, styles, plots, labels, ofh)

    colors = ["black", "blue", "green"]
    styles = ["--", "--", "--"]
    plots, labels, maxSens2 = plotRoc("CROP-IT", guideValidOts, cropitPredOts, colors, styles, plots, labels, ofh, isCropit=True)

    plt.legend(plots,
           labels,
           loc='lower right',
           ncol=1,
           fontsize=12)

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
    print "wrote data to %s" % ofh.name

main()
