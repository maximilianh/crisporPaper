import os, logging
logging.basicConfig(loglevel=logging.INFO)
from os.path import isfile, splitext, join
from annotateOffs import *
from collections import defaultdict

logging.basicConfig(loglevel=logging.INFO)

from scipy.stats import linregress, pearsonr, spearmanr

import matplotlib
matplotlib.use('Agg')
#matplotlib.rcParams['pdf.fonttype'] = 42
import matplotlib.pyplot as plt
import numpy as np

# normalize all scores and KO activity values to percent-rank ?
NORMALIZE = False

scoreCorrFh = None

scoreTypes = ["svm", "doench", "ssc", "chariRaw", "finalGc6", "finalGc2"]
scoreDescs = {
    "svm" : "Wang Score",
    "doench" : "Doench Score",
    "ssc" : "Xu Score",
    "chariRaw" : "Chari Score",
    "finalGc6" : "Ren, last 6 bp GC",
    "finalGc2" : "Farboud, last 2 bp GC"
}

datasetDescs = {
    "xu2015Train": "Wang/Xu Training",
    "doench2014-Hs": "Doench Human",
    "doench2014-Mm": "Doench Mouse", 
    "doench2014-CD33Exon2": "Doench CDS33 Exon2", 
    "varshney2015": "Varshney Zebra",
    "gagnon2014": "Gagnon Zebra",
    "xu2015": "Xu Validation",
    "ren2015": "Ren Drosophila",
    "farboud2015": "Farboud C. elegans",
    "chari2015Train": "Chari Human"
}

def plotScores(ax, scores, guideFreqs, scoreType, extFname, annotate, diam, isXu, scoreCorrFh, doLegend=False):
    " create scatter plot "
    data = defaultdict(list)
    names = []
    regrX = []
    regrY = []
    #print "REMOVING points where mod frequency is 0"
    dataset = basename(extFname).split(".")[0]

    for extSeq, (guideName, modFreq) in guideFreqs.iteritems():
        y = modFreq
        #if y==0.0:
            #continue
        x = scores[extSeq][scoreType]

        if annotate:
            names.append( (x, y, guideName) )

        if not isXu or not guideName.startswith("AAVS"):
            regrX.append(x)
            regrY.append(y)

        # just for plot: adding jitter for a scoretype with many identical scores
        if scoreType.startswith('finalGc'):
            x -= random.random()*0.5

        if isXu:
            if guideName.startswith("AR"):
                data["AR"].append( (x, y) )
            elif guideName.startswith("FOX"):
                data["FOXA1"].append( (x, y) )
            elif guideName.startswith("AAVS"):
                data["AAVS"].append( (x, y) )
            else:
                assert(False)
        else:
            data["data"].append( (x,y) )

    # only used during testing
    if NORMALIZE:
        newData  = {}
        for dataName, tuples in data.iteritems():
            xVals, yVals = zip(*tuples)
            xVals = convToRankPerc(xVals)
            yVals = convToRankPerc(yVals)
            dataPoints = zip(xVals, yVals)
            newData[dataName] = dataPoints
        data = newData
    # END - only used during testing

    if isXu:
        markers = {"AR" : "o", "FOXA1": "^", "AAVS":"x"}
    else:
        markers = {"data" : "o"}

    figs = []
    labels = []
    for title, tuples in data.iteritems():
        fig = ax.scatter(*zip(*tuples), \
            alpha=.5, \
            marker=markers[title], \
            linewidth=0, \
            s=diam)
        figs.append(fig)
        labels.append(title)

    if len(data)!=1:
        ax.legend(figs,
               labels,
               scatterpoints=1,
               loc='upper left',
               ncol=1,
               fontsize=10)

    if not NORMALIZE:
        slope, intercept, r_value, p_value, std_err = linregress(regrX,regrY)
        print "score type %s: Pearson R %f" % (scoreType, r_value)
        line = slope*np.asarray(regrX)+intercept
        ax.plot(regrX,line, linestyle='-', color="orange")
        #ax.annotate("R=%0.3f (p=%0.4f)" % (r_value, p_value), xy=(0.50,0.10), fontsize=8, xycoords='axes fraction')

        pearR, pearP = pearsonr(regrX, regrY)
        spearR, spearP = spearmanr(regrX, regrY)
        ret = pearR

        #row  = [scoreType, dataset, "%0.3f" % pearR, "%0.3f" % spearR]
        #scoreCorrFh.write("\t".join(row)+"\n")

        ax.annotate(r'Pearson R = %0.3f (p %0.3f)' % (pearR, pearP) + '\n' + r'Spearman $\rho$ = %0.3f (p %0.3f)' % (spearR, spearP), xy=(0.40,0.05), fontsize=9, xycoords='axes fraction')

    if annotate:
        for x, y, label in names:
           ax.annotate(
              label, fontsize=7, rotation=30, ha="right", rotation_mode="anchor",
              xy = (x, y), xytext = (0,0), alpha=0.9,
              textcoords = 'offset points', va = 'bottom')
    return ret


def extendTabAddScores(extFname, scores, scoreNames, outFname):
    " add columns for efficiency scores and write to extFname "
    #outFname = splitext(extFname)[0]+".scores.tab"
    ofh = None
    for row in iterTsvRows(extFname):
        if ofh==None:
            ofh = open(outFname, "w")
            ofh.write("\t".join(row._fields))
            ofh.write("\t")
            ofh.write("\t".join(scoreNames))
            ofh.write("\n")

        ofh.write("\t".join(row)+"\t")

        rowScores = []
        for name in scoreNames:
            rowScores.append(scores[row.extSeq][name])
        ofh.write("\t".join([str(x) for x in rowScores]))
        ofh.write("\n")
    ofh.close()
    print "wrote data to %s" % ofh.name

def plotDataset(datasetName, ax, db, title, yLabel="Knock-out efficiency", isXu=False, annotate=False, \
        diam=30, addColLabels=False, ylim=None):
    print ("plotting %s" % datasetName)
    scores, freqs = parseEffScores(datasetName, db)
    # create tsv
    extFname = join("effData/"+datasetName+".ext.tab")
    extendTabAddScores(extFname, scores, ["ssc", "doench", "svm", "chariRaw"], "out/%s-compEffData.tsv" % datasetName)

    ax[0].set_ylabel(yLabel)

    global scoreCorrFh

    corrs = []
    for index, scoreType in enumerate(scoreTypes):
        doLegend = (index==0)
        pearR = plotScores(ax[index], scores, freqs, scoreType, extFname, annotate, diam, isXu, scoreCorrFh, doLegend=doLegend)
        corrs.append(pearR)

    row = [datasetDescs[datasetName]]
    row.extend(corrs)
    row = [str(x) for x in row]
    scoreCorrFh.write("\t".join(row)+"\n")

    ax[0].set_xlim(0, 1.0)
    ax[1].set_xlim(0, 1.0)

    if ylim is not None:
        for axObj in ax:
            axObj.set_ylim(*ylim)

    if addColLabels:
        ax[0].set_title("SVM score from Wang et al. 2014")
        ax[1].set_title("Regr Score from Doench et al. 2014")
        ax[2].set_title("Regr Score from Xu et al. 2015")
        ax[3].set_title("SVM Score from Chari et al. 2015")
        ax[4].set_title("last 6bp GC count, Ren 2015 (+/-0.25")
        ax[5].set_title("last 2bp GC count, Farboud 2015 (+/-0.25")

    # put the row desc into the left border
    # http://stackoverflow.com/questions/25812255/row-and-column-headers-in-matplotlibs-subplots
    ax[0].annotate(title, xy=(0, 0.5), xytext=(-ax[0].yaxis.labelpad - 5, 0), \
       xycoords=ax[0].yaxis.label, textcoords='offset points', \
       #textcoords='offset points', \
       size='medium', ha='right', va='top')

def plotLargeScale():
    # large-scale studies to train the scoring models
    plotFname = "out/compEffScores-train.pdf" 

    fig, axArr = plt.subplots(4, 6, sharey="row")
    fig.set_size_inches(35,20)
    plotDataset("xu2015Train", axArr[0], "hg19", "Wang 2014\nhuman, HL-60 & KBM\nAs used by Xu 2015", diam=2, addColLabels=True, yLabel="rel. sgRNA abundance")
    plotDataset("doench2014-Hs", axArr[1], "hg19", "Doench 2014\nhuman MOLM13, NB4, TF1", diam=2)
    plotDataset("doench2014-Mm", axArr[2], "mm9", "Doench 2014\nmouse EL4", diam=2)
    plotDataset("chari2015Train", axArr[3], "hg19", "Chari 2015\nhuman K562", diam=2, yLabel="Mutation Rate", ylim=(0,2.0))

    fig.tight_layout()
    fig.subplots_adjust(left=0.15, top=0.95)
    fig.savefig(plotFname, format = 'pdf')
    fig.savefig(plotFname.replace(".pdf", ".png"))
    print "wrote plot to %s, added .png" % plotFname
    plt.close()

def main():
    #"XuData/modFreq.tab" 
    #extendTabAddContext("temp.tab")
    global scoreCorrFh
    scoreCorrFh = open("out/effScoreComp.tsv", "w")
    scoreCorrFh.write("dataset\t%s\n" % "\t".join([scoreDescs[st] for st in scoreTypes]))


    plotLargeScale()
    # small-scale studies
    fig, axArr = plt.subplots(11, 6, sharey="row")
    fig.set_size_inches(35,11*4)

    plotFname = "out/compEffScores-valid.pdf"

    plotDataset("xu2015", axArr[0], "hg19", "Xu 2015 validation\nFOX/AR: Western Blot\nAAVS1: T7 assay\nRegression/Correlations only for FOX/AR", isXu=True, addColLabels=True)
    plotDataset("doench2014-CD33Exon2", axArr[1], "hg19", "Doench\nCDS33 Exon2")
    plotDataset("varshney2015", axArr[2], "danRer10", "Varshney 2015\nZebrafish egg RNA injection")
    plotDataset("gagnon2014", axArr[3], "danRer10", "Gagnon 2014\nZebrafish egg RNA injection")
    plotDataset("ren2015", axArr[4], "dm3", "Ren 2015\nDrosophila injection")
    plotDataset("farboud2015", axArr[5], "ce6", "Farboud 2015\nC. elegans injection")

    axStart = 6
    chariCells = ["293T", "A549", "HepG2", "K562", "PGP1iPS", "SKNAS", "U2OS"]
    for i in range(0, len(chariCells)):
        cell = chariCells[i]
        plotDataset("chari2015Valid_"+cell, axArr[axStart+i], "hg19", "Chari 2015\nhuman, "+cell)

    fig.tight_layout()
    fig.subplots_adjust(left=0.15, top=0.95)
    fig.savefig(plotFname, format = 'pdf')
    fig.savefig(plotFname.replace(".pdf", ".png"))
    print "wrote plot to %s, added .png" % plotFname
    plt.close()


    print "wrote score R summary to %s" % scoreCorrFh.name

main()
