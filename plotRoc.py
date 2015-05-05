import glob, sys
from annotateOffs import *

import numpy as np
import matplotlib.pyplot as plt
from os.path import basename, splitext

#minFrac= 0.02
#minFrac = 0.01

onlyAlt = False

if len(sys.argv)>1:
    altPamCutoff = float(sys.argv[1])
else:
    altPamCutoff = None

def parse(fname):
    """ parse the annotated validated off-target table and return as dict
    guideSeq -> otSeq -> otScore. Also guideName -> guideSeq.
    """
    otScores = defaultdict(dict)
    for row in iterTsvRows(fname):
        if not row.otSeq[-2:] in ["AG", "GA"]:
            if onlyAlt:
                continue
        otScores[row.guideSeq][row.otSeq] = float(row.readFraction)
    return otScores

def parseCrispor(dirName):
    " parse rispor output files, return as dict guideSeq -> ot seq -> otScore "
    predScores = defaultdict(dict)
    targetSeqs = {}
    for fname in glob.glob(dirName+"/*.tsv"):
        for row in iterTsvRows(fname):
            guideName = splitext(basename(fname))[0]
            predScores[row.guideSeq][row.offtargetSeq] = float(row.offtargetScore)
            targetSeqs[guideName] = row.guideSeq
    return predScores, targetSeqs

def parseMit(dirName, targetSeqs):
    " parse the MIT csv files, return a dict with guideSeq -> otSeq -> otScore "
    fnames= glob.glob(dirName+"/*.csv")
    data = defaultdict(dict)
    for fname in fnames:
        guideName = fname.split("/")[1].split(".")[0]
        study = guideName.split("_")[0]
        #if study in ignoreStudies:
            #continue
        if guideName not in targetSeqs:
            print "MIT off-target data without bench data: %s" % guideName
            continue
        for line in open(fname):
            if line.startswith("guide"):
                continue
            fs = line.split(", ")
            otSeq = fs[4]
            score = fs[6]
            if fs[8]=="True":
                # ontarget
                continue
            guideSeq = targetSeqs[guideName]
            data[guideSeq][otSeq]=float(score)
    return data


def getRocValues(toolName, guideValidOts, guidePredOts, minReadFrac):
    " return a list of (sens, fdr) tuples for a ROC curve plotting "
    # keep only the validated off-targets with read fraction > minCutoff
    validOts = set()
    for guideSeq, validOtSeqs in guideValidOts.iteritems():
        for seq, readFrac in validOtSeqs.iteritems():
            if readFrac > minReadFrac:
                validOts.add(seq)

    sensList = []
    fdrList = []

    for cutoff in [0.001, 0.002, 0.003, 0.005, 0.007, 0.009, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06, 0.07, 0.08, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, 3, 4, 5, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30, 35, 40, 50, 60]:
        predOts = set()
        allOts = set()
        allOts.update(validOts) # make a copy of the elements

        for guideSeq, predSeqScores in guidePredOts.iteritems():
            # get the predicted sequences over the off-target score cutoff
            for predSeq, seqScore in predSeqScores.iteritems():
                allOts.add(predSeq)

                # check if alternative PAM
                if predSeq[-2:] in ["AG", "GA"]:
                    if altPamCutoff!=None and seqScore < altPamCutoff:
                        continue
                else:
                    if onlyAlt:
                        continue

                if seqScore > float(cutoff):
                    predOts.add(predSeq)

        notPredOts = allOts - predOts
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
        print "\t".join(row)
        sys.stdout.flush()

    return sensList, fdrList, validOts

def plotRoc(prefix, guideValidOts, guidePredOts, colors, styles, plots, labels):
    i= 0
    maxSens = 0
    for minFrac in [0.0, 0.001, 0.01]:
        sensList, fdrList, validSeqs = getRocValues(prefix, guideValidOts, guidePredOts, minFrac)
        plotLabel = prefix+", cleavage > %0.1f%%, %d off-targets" % ((minFrac*100), len(validSeqs))
        p, = plt.plot(fdrList, sensList, ls=styles[i], color=colors[i]) # NB: comma!
        plots.append(p)
        labels.append(plotLabel)
        maxSens = max(maxSens, max(sensList))
        i+=1
    return plots, labels, maxSens

guideValidOts = parse("annotFiltOfftargets.tsv")
guidePredOts, targetSeqs = parseCrispor("crisporOfftargets")
mitPredOts = parseMit("mitOfftargets", targetSeqs)

headers = ["readFrac", "cutoff", "sens", "spec", "TP", "FP", "FN", "TN"]
print "\t".join(headers)

plots = []
labels = []

plots = []
labels = []

styles = ["-", "-", "-"]
colors = ["black", "blue", "green"]
plots, labels, maxSens1 = plotRoc("Crispor", guideValidOts, guidePredOts, colors, styles, plots, labels)

colors = ["black", "blue", "green"]
styles = [":", ":", ":"]
plots, labels, maxSens2 = plotRoc("MIT", guideValidOts, mitPredOts, colors, styles, plots, labels)

plt.legend(plots,
       labels,
       loc='lower right',
       ncol=1,
       fontsize=12)

plt.xlabel("False positive rate")
plt.ylabel("True positive rate")
ax = plt.gca()
ax.axhline(y=maxSens1, ls=":", color="k")
plt.text(0, maxSens1, "max = %0.2f"%maxSens1)
plt.ylim(0,1.0)
plt.xlim(0,1.0)

outfname = "roc.pdf"
plt.savefig(outfname)
print "wrote %s" % outfname

outfname = "roc.png"
plt.savefig(outfname)
print "wrote %s" % outfname
