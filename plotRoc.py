# plot the ROC curve of the off-targets when varying the off-target score
import glob, sys, logging
from annotateOffs import *

import numpy as np
import matplotlib.pyplot as plt
from os.path import basename, splitext, isfile
from sklearn.metrics import roc_curve, auc, precision_recall_curve

import crisporOtScores

# a supplemental version of this figure can be created by specifying "supp" as the 
# first argument. The supp version has maxMismatches increased to 6 and adds an 
# additional dataset that includes the two outliers that are usually removed
# from the analysis by filtAnnotateOfftargets.py
makeSupp = False

# we only look at off-targets with a certain number of mismatches
# otherwise the ROC curve would not go up to 1.0 as MIT can only search for 4 MMs
maxMismatches = 4

# only look at alternative PAMs, can be used to determine best cutoff for the alternative PAMs
#onlyAlt = True
onlyAlt = False

# for the ROC curve, we only analyze off-targets with certain PAM sites 
# assuming that no software can find the relatively rare PAM sites
# that are not GG/GA/AG
validPams = ["GG", "GA", "AG"]
if onlyAlt:
    validPams = ["GA", "AG"]


#if len(sys.argv)>1:
    #altPamCutoff = float(sys.argv[1])
#else:
# altPamCutoff = None

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

def plotRoc(prefix, expFreqs, predScores, color, style, plots, labels, allOts=None):
    """ plot ROC curve and write annotation to out file. Add result to plots and labels lists of objects.  """
    if allOts is None:
        allOts = predScores.keys()

    print "plotting ROC for %s" % prefix

    # use all off-targets, don't use the modification frequency to restrict the elements
    minFrac = 0.0

    yVals = []
    yTrue = []
    for otSeq in allOts:
        if "MIT" in prefix:
            # the MIT site gives us no score for many off-targets
            # so we're setting it to 0.0 for these
            # it's not entirely correct, as we should somehow treat these as "missing data"
            # this leads to a diagonal line in the plot... not sure how to avoid that
            otScore = predScores.get(otSeq, 0.0)
        else:
            otScore = predScores[otSeq]
        #print otScore
        assert(otScore is not None)
        yVals.append(otScore)
        # only ignore off-targets with a mod freq of 0.0 (can appear in Hsu/Cho data)
        if expFreqs.get(otSeq, 0.0) > minFrac:
            yTrue.append(1)
        else:
            yTrue.append(0)

    assert(len(yTrue)==len(yVals))
    assert(None not in yTrue)
    assert(None not in yVals)
    print "%s, minFreq %f, number of scored elements: %d "% (repr(prefix), minFrac, len(yVals))
    print("number of positives: %d" % (sum(yTrue)))
    #print "yTrue", yTrue
    #print "yVals", yVals
    fpr, tpr, thresholds = roc_curve(yTrue, yVals)
    #fpr, tpr, thresholds = precision_recall_curve(yTrue, yVals)
    scoreName = prefix.split()[0]
    if scoreName == "CFD":
        ofh = open("out/%s-mm%d-%0.3f-rocThresholds.tsv" % (scoreName, maxMismatches, minFrac), "w")
        ofh.write("fpr\ttpr\tthreshold\n")
        for f, t, score in zip(fpr, tpr, thresholds):
            row = [f, t, score]
            row = [str(x) for x in row]
            ofh.write("\t".join(row)+"\n")
        ofh.close()
        print "Wrote ROC thresholds to %s" % ofh.name
    roc_auc = auc(fpr, tpr)
    print "AUC=%f" % roc_auc
    plotLabel = prefix + ", AUC %.2f" % roc_auc
    #plotLabel = prefix
    # plotLabel = prefix+", mod. freq. > %0.1f%%: AUC %0.2f" % ((minFrac*100), roc_auc)
    p, = plt.plot(fpr, tpr, ls=style, color=color) # NB: keep the comma!
    plots.append(p)
    labels.append(plotLabel)
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
    ofh.write("\t".join(["guide", "offt", "modFreq", "MIT", "hsu", "cfd"])+"\n")
    for guideSeq, otFreqs in guideValidOts.iteritems():
        for otSeq, otFreq in otFreqs.iteritems():
            mitScore = crisporOtScores.calcMitScore(guideSeq, otSeq)
            hsuScore = crisporOtScores.calcHsuSuppScore2(guideSeq, otSeq)
            cfdScore = crisporOtScores.calcCfdScore(guideSeq, otSeq)
            row = [guideSeq, otSeq, otFreq, mitScore, hsuScore, cfdScore]
            row = [str(x) for x in row]
            ofh.write("\t".join(row)+"\n")
    print "wrote %s" % ofh.name

def parseOfftargets(fname, maxMismatches, onlyAlt, validPams):
    """ parse the annotated validated off-target table and return as dict
    guideSeq -> otSeq -> modifFreq and another dict guideName -> guideSeq
    """
    otScores = defaultdict(dict)
    guideSeqs = dict()
    print "parsing %s" % fname
    skipCount = 0
    for row in iterTsvRows(fname):
        #print fname, int(row.mismatches), maxMismatches
        if int(row.mismatches)>maxMismatches:
            print "skip", row
            skipCount += 1
            continue
        if validPams!=None and not row.otSeq[-2:] in validPams:
            print "not using off-target %s/%s, PAM is not in %s" % (row.name, row.otSeq, validPams)
            continue

        guideSeqs[row.name] = row.guideSeq
        if onlyAlt and not row.otSeq[-2:] in ["AG", "GA"]:
            continue
        otScores[row.guideSeq][row.otSeq] = float(row.readFraction)
    print "Skipped %d rows with more than %d mismatches" % (skipCount, maxMismatches)
    return otScores, guideSeqs

def main():

    if len(sys.argv)>1:
        assert(sys.argv[1]=="supp")
        if sys.argv[1]=="supp":
            global makeSupp
            makeSupp = True

    #if makeSupp:
        #global maxMismatches
        #maxMismatches = 6
        #print "Setting max mismatches to 6"

    guideValidOts, guideSeqs = parseOfftargets("out/annotFiltOfftargets.tsv", maxMismatches, onlyAlt, validPams)

    mitPredOts = parseMit("mitOfftargets", guideSeqs)

    tmpFname = "/tmp/crisporOffs-%d.pickle" % maxMismatches
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

    crisporScores = collapseDicts(crisporPredOts)
    print "Loaded %d CRISPOR predicted off-targets" % len(crisporScores)

    validOffts = collapseDicts(guideValidOts)

    print "calculating CFD scores"
    cfdPredOts = calcOtScores(crisporPredOts, calcCfdScore)
    cfdScores = collapseDicts(cfdPredOts)
    plots, labels = plotRoc("CFD score", validOffts, cfdScores, "orange", "-", plots, labels)

    if makeSupp:
        # do not remove the two Tsai et al outliers
        guideValidOtsFull, _ = parseOfftargets("out/annotOfftargets.tsv", maxMismatches, onlyAlt, validPams)
        validOfftsFull = collapseDicts(guideValidOtsFull)
        plots, labels = plotRoc("CFD, with two outliers", validOfftsFull, cfdScores, "grey", ":", plots, labels)
        #print len(validOfftsFull)
        # assert ("AACCACCACTCCAGACTAAAGAG" in validOfftsFull)
        #print "XX", len(validOfftsFull)

    print "plotting other ROCs"
    plots, labels = plotRoc("MIT score", validOffts, crisporScores, "darkblue", "-", plots, labels)
    if not onlyAlt:
        plt.annotate('FPR=0.43, TPR=0.98:\nCFD score = 0.023', xy=(0.44, 0.985), xytext=(0.5, 1.04),
                arrowprops=dict(facecolor='black', arrowstyle="->"), annotation_clip=False)

    mitScores = collapseDicts(mitPredOts)
    plots, labels = plotRoc("MIT Website", validOffts, mitScores, "blue", "-", plots, labels, allOts=crisporScores.keys())

    cropitScores = collapseDicts(calcOtScores(crisporPredOts, calcCropitScore))
    plots, labels = plotRoc("Cropit score", validOffts, cropitScores, "red", "-.", plots, labels)

    ccTopPredOts = calcOtScores(crisporPredOts, calcCcTopScore)
    ccTopScores = collapseDicts(ccTopPredOts)
    plots, labels = plotRoc("CCTop score", validOffts, ccTopScores, "green", ":", plots, labels)

    #hsuScores = collapseDicts(calcOtScores(crisporPredOts, crisporOtScores.calcHsuSuppScore2))
    #plots, labels = plotRoc("Hsu score", validOffts, hsuScores, "grey", "--", plots, labels)

    #plots, labels = plotRoc("CROP-IT", crisporScores, validOffts, cropitPredOts, colors, styles, plots, labels, isCropit=True)

    plt.legend(plots,
           labels,
           loc='lower right',
           ncol=1,
           fontsize=13)

    plt.xlabel("False positive rate (FPR)")
    plt.ylabel("True positive rate (TPR)")

    plt.ylim(0,1.0)
    plt.xlim(0,1.0)

    outfname = "out/roc.pdf"
    if makeSupp:
        outfname = "out/roc-supp.pdf"

    plt.savefig(outfname)
    print "wrote %s" % outfname

    outfname = outfname.replace(".pdf", ".png")
    plt.savefig(outfname)
    print "wrote %s" % outfname
    #print "wrote data to %s" % ofh.name

    if not makeSupp:
        compHsuMit(guideValidOts, guideSeqs, "out/hsuVsMit.tsv")


main()
