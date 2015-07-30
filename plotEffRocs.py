import os, logging
from os.path import isfile, splitext, join, basename
from annotateOffs import *
from collections import defaultdict
from sklearn.metrics import roc_curve, auc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

def main():
    #"XuData/modFreq.tab" 
    #extendTabAddContext("temp.tab")
    fig, axArr = plt.subplots(7, 1)
    fig.set_size_inches(5, 30)

    #fnames = glob.glob("effData/*.ext.tab")
    datasets = ["xu2015Train", "doench2014-Hs", "doench2014-Mm", "varshney2015", "gagnon2014", "xu2015", "ren2015"]
    for row, dataset in enumerate(datasets):
        #dataset = basename(fname).split(".")[0]
        seqScores, guideFreqs = parseEffScores(dataset)
        ax = axArr[row]
        ax.set_title(dataset)
        ax.set_xlabel('False Positive Rate')
        ax.set_ylabel('True Positive Rate')
        ax.set_xlim([0.0, 1.0])
        ax.set_ylim([0.0, 1.05])

        # get minimum eff score
        effVals = [koScore for name, koScore in guideFreqs.values()]
        effCutoff = np.percentile(effVals, 80)

        scoreTypes = ["svm", "doench", "ssc", "chariRaw", "finalGc6"]
        colors = ["red", "black", "blue", "green", "orange"]
        lineTypes = ["-", "--", "-.", ":", "-", "--"]

        for scoreIdx, scoreType in enumerate(scoreTypes):
            xVals = []
            classNames = []
            for seq, scores in seqScores.iteritems():
                guideName, modFreq = guideFreqs[seq]
                xVals.append(scores[scoreType])

                if modFreq > effCutoff:
                    classNames.append(1)
                else:
                    classNames.append(0)

            fpr, tpr, thresholds = roc_curve(classNames, xVals)
            roc_auc = auc(fpr, tpr)
            print scoreType, roc_auc

            ax.plot(fpr, tpr, label='{0} (AUC = {1:0.2f})'.format(scoreType, roc_auc), \
                color=colors[scoreIdx], linestyle=lineTypes[scoreIdx])
        ax.legend(loc="lower right", fontsize=9)

    plotFname = "out/effRoc.pdf"

    plt.tight_layout()
    plt.savefig(plotFname)
    plt.savefig(plotFname.replace(".pdf", ".png"))
    print "wrote plot to %s and .png" % plotFname
main()
