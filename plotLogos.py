# plot weblogos of the regression coefficients in the doench,chari and xu/wang
# datasets

import os
import numpy as np
from sklearn.linear_model import *

from annotateOffs import *
import plotSeqLogo

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def parseTab(fname, perc):
    " return all seqs with scores that are higher than the perc percentile "
    xList = []
    yList = []
    seqs = []
    scores = []
    for row in iterTsvRows(open(fname)):
        seq = row.seq[:20]
        seqs.append(seq)
        score = float(row.modFreq)
        scores.append(score)

    #sortScores = list(sorted(scores))
    #cutoff = sortScores[-100]
    #cutoff = np.percentile(scores, perc)
    #print "cutoff", cutoff

    #filtSeqs = []
    #for seq, score in zip(seqs, scores):
        #if score > cutoff:
            #filtSeqs.append(seq)
    #print "sequences", len(filtSeqs)
    return seqs, scores

#def makeLogo(seqs, title, fname):
    #ofh = open("/tmp/seqs.fa", "w")
    #for i, seq in enumerate(seqs):
        #ofh.write(">%d\n%s\n" % (i, seq))
    #ofh.close()
    #print ofh.name
#
    #cmd = "weblogo -t '%s' < /tmp/seqs.fa -F png -U probability > %s" % (title, fname)
    #os.system(cmd)
    #print "wrote %s" % fname
#
    #os.remove("/tmp/seqs.fa")

def main():
    datasets = ["doench2014-Hs", "xu2015Train", "chari2015Train", "varshney2015", "gagnon2014"]
    #datasets = ["doench2014-Hs", "xu2015Train"]
    fnames = []
    fig, axArr = plt.subplots(len(datasets), sharex=True)
    fig.set_size_inches(10,10)
    for plotIdx, dataset in enumerate(datasets):
        ax = axArr[plotIdx]
        ax.set_xlim(0,80)
        ax.set_title(dataset)
        seqs, scores = parseTab("effData/%s.tab" % dataset, 75)
        vecs = seqsToVecs(seqs)
        #regr = LinearRegression(normalize=True)
        regr = Lasso(alpha=0.001)
        #regr = ElasticNet(alpha=0.0)
        regr.fit(vecs, scores)
        #fname = "out/logo-%s.png" % dataset
        #fnames.append(fname)
        #ddmakeLogo(seqs, dataset, fname)
        minCoef, maxCoef = min(regr.coef_), max(regr.coef_)
        print "before"
        printCoef(regr.coef_)
        normCoef = [x*(1.0/maxCoef) for x in regr.coef_]
        print "after norm:"
        printCoef(normCoef)
        logoDicts = vecToSeqDicts(normCoef)
        plotSeqLogo.plotLogo(fig, ax, logoDicts)

    plt.savefig("out/logos.png")
    plt.close()

    #cmd = "convert -append %s out/logos.png" % " ".join(fnames)
    #os.system(cmd)


main()
