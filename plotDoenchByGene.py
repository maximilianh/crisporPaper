import os, logging
from annotateOffs import *
from collections import defaultdict
#from scipy.stats import linregress

from scipy.stats import linregress

import matplotlib.pyplot as plt
import numpy as np

def parseDoenchByGene(fname):
    " return a dict geneName -> seq -> activity "
    data = defaultdict(dict)
    for row in iterTsvRows(fname):
       gene = row.guide.split("-")[0]
       data[gene][row.extSeq] = float(row.modFreq)
    return data

def main():
    #plt.figure(figsize=(8,10))
    fig, axArr = plt.subplots(6, 4, sharex="col", sharey="row")
    fig.set_size_inches(20,15)

    actData = parseDoenchByGene("effData/doench2014-Mm.ext.tab")
    for plotRow, gene in enumerate(actData):
        print plotRow, gene
        seqKoFreqs = actData[gene]
        scores = calcEffScores(seqKoFreqs)

        for plotCol, scoreType in enumerate(["svm", "doench", "ssc", "chariRaw"]):
            xVals = []
            yVals = []
            for extSeq, koScore in seqKoFreqs.iteritems():
                seqScores = scores[extSeq]
                xVals.append(seqScores[scoreType])
                yVals.append(koScore)

            ax = axArr[plotRow][plotCol]
            ax.scatter(xVals, yVals, alpha=0.5, linewidths=0)
            ax.set_xlabel("%s score" % scoreType)
            ax.set_ylabel("Knock-out assay result")
            ax.set_ylim(0,1.0)
            if plotCol==0:
                ax.set_title(gene)

            slope, intercept, r_value, p_value, std_err = linregress(xVals,yVals)
            print 'r^2 value', r_value**2
            print 'p_value', p_value
            print 'standard deviation', std_err
            line = slope*np.asarray(xVals)+intercept
            ax.plot(xVals,line, linestyle='--', color="orange")
            ax.annotate("r^2=%0.3f\npVa=%0.4f" % (r_value**2, p_value), xy=(0.75,0.05), fontsize=8, xycoords='axes fraction')

    outFname = "out/doenchByGene.pdf"
    fig.tight_layout()
    plt.savefig(outFname)
    plt.savefig(outFname.replace(".pdf", ".png"))
    print "wrote plot to %s and .png" % outFname

main()
