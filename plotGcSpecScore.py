#scatter plot of GC content versus specificity score
import os, logging
from annotateOffs import *
from collections import defaultdict
import matplotlib.pyplot as plt
import string

def parseGcSpec(fname):
    """ return dicts with guideName -> gcConts and guideName -> specScores 
    and guideName -> labels
    """
    gcConts = defaultdict(list)
    specScores = defaultdict(list)
    labels = defaultdict(list)
    doneGuides = set()

    for row in iterTsvRows(fname):
        study, guideName = string.split(row.name, "_", maxsplit=1)
        #guideGc = (row.seq[:20].count("G") + row.seq[:20].count("C")) / 20.0
        if row.guideSeq in doneGuides:
            continue
        doneGuides.add(row.guideSeq)

        gcConts[study].append(gcCont(row.guideSeq[:20]))
        specScores[study].append(int(row.guideSpecScore4MM))
        labels[study].append(guideName)
    return (gcConts, specScores, labels)

def main():
    gcConts, specScores, labels = parseGcSpec("out/annotOfftargets.tsv")

    studies = sorted(gcConts.keys())

    colors  = ["green", "blue", "green", "blue", "red", "grey", "orange", "blue", "orange", "red", "magenta", "blue"]
    markers = ["o", "s", "+", ">", "<", "o", ".", "o", "x", "+", ".", "<"]
    for i, study in enumerate(studies):
        xVals = gcConts[study]
        yVals = specScores[study]
        studyFig = plt.scatter(xVals, yVals, \
            alpha=0.9, \
            marker=markers[i], \
            s=30, \
            color=colors[i])
        for x, y, label in zip(xVals, yVals, labels[study]):
            plt.annotate( \
              label, fontsize=9, rotation=0, ha="right", rotation_mode="anchor", \
              xy = (x, y), xytext = (-2,-1), alpha=1.0, textcoords = 'offset points', va = 'bottom')
    plt.xlabel("GC content of guide sequence")
    plt.ylabel("Specificity score of guide sequence")

    outfname = "out/gcSpecScores.pdf"
    plt.savefig(outfname, format = 'pdf')
    plt.savefig(outfname.replace(".pdf", ".png"))
    print "wrote %s" % outfname

main()
