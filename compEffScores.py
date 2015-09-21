import os, logging
logging.basicConfig(loglevel=logging.INFO)
from os.path import isfile, splitext, join
from annotateOffs import *
from collections import defaultdict

logging.basicConfig(loglevel=logging.INFO)

from scipy.stats import linregress, pearsonr, spearmanr, mannwhitneyu, rankdata

import matplotlib
matplotlib.use('Agg')
#matplotlib.rcParams['pdf.fonttype'] = 42
import matplotlib.pyplot as plt
import numpy as np

# normalize all scores and KO activity values to percent-rank ?
NORMALIZE = False

scoreCorrFh = None

# the last two score types are not written to the tab-sep file
scoreTypes = ["wangOrig", "doench", "ssc", "chariRaw", "crisprScan", "drsc", "fusi", "finalGc6", "finalGg"]

scoreDescs = {
    "wang" : "Wang Score2",
    "wangOrig" : "Wang Score",
    "doench" : "Doench Score",
    "ssc" : "Xu (Wang) Score",
    "chariRaw" : "Chari Score",
    "chariRank" : "Chari Rank Score",
    "chariRaw" : "Chari Score",
    "finalGc6" : "Ren: last 6 bp GC>4",
    #"finalGc2" : "Farboud-like, last 2 bp GC",
    "finalGg" : "Farboud: ends with GG",
    "myScore" : "Max Score",
    "crisprScan" : "Morenos-Matos Score",
    "fusi" : "Fusi/Doench Score",
    "drsc" : "Housden Score",
    "oof"      : "Bae Out-of-Frame Score"
}

# labels for the different score types
scoreLabels = {
        "wang" : "SVM score2 from Wang et al. 2014",
        "wangOrig" : "SVM score from Wang et al. 2014",
        "doench" : "Regr Score from Doench et al. 2014",
        "chariRank" : "SVM Rank Score Rank from Chari et al. 2015",
        "chariRaw" : "SVM Score from Chari et al. 2015",
        "ssc" : "Regr Score from Xu et al. 2015",
        "fusi" : "Regr Score from Fusi/Doench et al. 2015",
        "drsc" : "Score from Housden et al. 2015",
        "crisprScan" : "Regr Score from Morenos-Matos et al. 2015",
        "finalGc6" : "last 6bp GC>4, Ren 2015 +/-0.25",
        "finalGg" : "last 2bp=GG , Farboud 2015 +/-0.25"
        }

datasetDescs = {
    "xu2015Train": "Wang/Xu",
    "xu2015TrainHl60": "Wang/Xu HL60",
    "xu2015TrainKbm7": "Wang/Xu KBM7",
    "doench2014-Hs": "Doench Human MOLM13/NB4/TF1",
    "doench2014-Mm": "Doench Mouse",
    "housden2015": "Housden Dros-S2R+",
    "doench2014-CD33Exon2": "Doench CD33 Exon2",
    "doench2014-CD33Exon3": "Doench CD33 Exon3",
    "doench2014-CD13Exon10": "Doench CD13 Exon10",
    "varshney2015": "Varshney Zebrafish",
    "gagnon2014": "Gagnon Zebrafish",
    "xu2015": "Xu Validation",
    "xu2015FOX-AR": "Xu KO Validation FOX/AR",
    "xu2015AAVS1": "Xu T7 Validation AAVS1",
    "ren2015": "Ren Drosophila",
    "farboud2015": "Farboud C. elegans",
    "morenoMateos2015": "Moreno-Mateos Zebrafish",
    "museumT7": "Concordet Indel T7",
    "museumIC50": "Concordet Indel IC50",
    "schoenig": "Schoenig Indel LacZ Rank",
    "eschstruth": "Eschstruth Indel T7 Rank",
    "chari2015Train": "Chari Human",
    "chari2015Train293T": "Chari 293T",
    "concordet2-Hs": "Conc2 Hs",
    "concordet2-Mm": "Conc2 Mm",
    "concordet2-Rn": "Conc2 Rn",
    "concordet2": "Concordet2 Hs/Mm/Rn",
    "chari2015TrainK562": "Chari Human K562"
}

def plotScores(ax, scores, guideFreqs, scoreType, annotate, diam, doLegend=False):
    " create scatter plot "
    regrX = []
    regrY = []
    plotX = []
    plotY = []

    for extSeq, (guideName, modFreq) in guideFreqs.iteritems():
        y = modFreq
        x = scores[extSeq][scoreType]

        regrX.append(x)
        regrY.append(y)

        # just for plot: adding jitter for a scoretype with many identical scores
        if scoreType.startswith('final'):
            x -= random.random()*0.25

        plotX.append(x)
        plotY.append(y)

    ax.scatter(plotX, plotY, alpha=.5, marker="o", s=diam, linewidth=0)

    if scoreType in ["wang", "wangOrig", "doench"]:
        ax.set_xlim(0, 1.0)
    elif scoreType=="chariRank":
        ax.set_xlim(0, 100.0)

    slope, intercept, r_value, p_value, std_err = linregress(regrX,regrY)
    print "score type %s: Pearson R %f" % (scoreType, r_value)
    line = slope*np.asarray(regrX)+intercept
    ax.plot(regrX,line, linestyle='-', color="orange")

    pearR, pearP = pearsonr(regrX, regrY)
    spearR, spearP = spearmanr(rankdata(regrX), rankdata(regrY))
    #mwU, mwP = mannwhitneyu(regrX, regrY)
    ret = spearR
    #ret = pearR

    #ax.annotate(r'Pearson R = %0.3f (p %0.3f)' % (pearR, pearP) + '\n' + r'Spearman $\rho$ = %0.3f (p %0.3f)' % (spearR, spearP) + "\nMann-Whitney U=%d (p=%0.3f)" % (int(mwU), mwP), xy=(0.40,0.08), fontsize=9, xycoords='axes fraction')
    ax.annotate(r'Pearson R = %0.3f (p %0.3f)' % (pearR, pearP) + '\n' + r'Spearman $\rho$ = %0.3f (p %0.3f)' % (spearR, spearP), xy=(0.40,0.08), fontsize=9, xycoords='axes fraction')

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

def plotDataset(datasetName, ax, title, yLabel="Knock-out efficiency", annotate=False, \
        diam=30, addColLabels=True, ylim=None, yTicks=None):
    global scoreCorrFh

    print ("plotting %s" % datasetName)
    scores, freqs = parseEffScores(datasetName)
    ax[0].set_ylabel(yLabel)

    corrs = []
    for index, scoreType in enumerate(scoreTypes):
        doLegend = (index==0)
        corr = plotScores(ax[index], scores, freqs, scoreType, annotate, diam, doLegend=doLegend)
        #if scoreType not in ["finalGc6", "finalGg"]:
        corrs.append(corr)

    # write to tab sep file for heatmap
    row = [datasetDescs[datasetName] + (" (%d)" % len(freqs)) ]
    row.extend(["%0.3f" % c for c in corrs])
    if scoreCorrFh is not None:
        scoreCorrFh.write("\t".join(row)+"\n")

    #if datasetName.startswith("doench"):
        #for a in ax:
            #a.set_ylim(0, 1.0)
    #ax[0].set_xlim(0, 1.0)
    #ax[1].set_xlim(0, 1.0)
    #ax[3].set_xlim(0, 100)

    if ylim is not None:
        for axObj in ax:
            axObj.set_ylim(*ylim)

    if addColLabels:
        for i in range(len(scoreTypes)):
            ax[i].set_title(scoreLabels[scoreTypes[i]])

    # put the row desc into the left border
    # http://stackoverflow.com/questions/25812255/row-and-column-headers-in-matplotlibs-subplots
    ax[0].annotate(title, xy=(0, 0.5), xytext=(-ax[0].yaxis.labelpad - 5, 0), \
       xycoords=ax[0].yaxis.label, textcoords='offset points', \
       #textcoords='offset points', \
       size='medium', ha='right', va='top')

    if yTicks:
        ax[0].set_yticks(yTicks)

def plotLargeScale(corrFname):
    # large-scale studies to train the scoring models
    global scoreCorrFh
    scoreCorrFh = open(corrFname, "w")
    scoreCorrFh.write("dataset\t%s\n" % "\t".join([scoreDescs[st] for st in scoreTypes]))

    plotFname = "out/compEffScores-train.pdf" 

    rowCount = 5
    fig, axArr = plt.subplots(rowCount, len(scoreTypes))
    fig.set_size_inches(len(scoreTypes)*5,rowCount*5)
    #plotDataset("xu2015Train", axArr[0], "hg19", "Wang 2014\nhuman HL-60\nAs used by Xu 2015", diam=2, addColLabels=True, yLabel="-log2 sgRNA fold change")
    plotDataset("xu2015TrainHl60", axArr[0], "Wang 2014\nhuman HL-60\nAs used by Xu 2015", diam=2, addColLabels=True, yLabel="-log2 sgRNA fold change")
    #plotDataset("xu2015TrainKbm7", axArr[1], "Wang 2014\nhuman KBM7\nAs used by Xu 2015", diam=2, addColLabels=True, yLabel="-log2 sgRNA fold change")
    plotDataset("doench2014-Hs", axArr[1], "Doench 2014\nhuman MOLM13, NB4, TF1", diam=2, yLabel="rank-percent", ylim=(0,1.0))
    plotDataset("doench2014-Mm", axArr[2], "Doench 2014\nmouse EL4", diam=2, yLabel="rank-percent", ylim=(0,1.0))
    plotDataset("chari2015Train293T", axArr[3], "Chari 2015\nhuman 293T", diam=2, yLabel="Mutation Rate", ylim=(0,2.0))
    plotDataset("chari2015TrainK562", axArr[4], "Chari 2015\nhuman K562", diam=2, yLabel="Mutation Rate")

    i = 5
    global datasetDescs
    #for dataset in ["doench2014-Hs-MOLM13_CD15","doench2014-Hs-NB4_CD13","doench2014-Hs-TF1_CD13", "doench2014-Hs-MOLM13_CD33","doench2014-Hs-NB4_CD33","doench2014-Hs-TF1_CD33"]:
    #for dataset in ["doench2014-Hs-MOLM13_CD15","doench2014-Hs-NB4_CD13","doench2014-Hs-TF1_CD13"]:
        #datasetDescs[dataset] = dataset
        #cellType = dataset.split("-")[-1].replace("_", " ")
        #plotDataset(dataset, axArr[i], "hg19", "Doench 2015\nhuman %s" % cellType, diam=2, yLabel="Fold Abundance")
        #i += 1

    fig.tight_layout()
    fig.subplots_adjust(left=0.15, top=0.95)
    fig.savefig(plotFname, format = 'pdf')
    fig.savefig(plotFname.replace(".pdf", ".png"))
    print "wrote plot to %s, added .png" % plotFname
    plt.close()

def plotSmallScale():
    # plot small-scale studies

    figCount = 15
    fig, axArr = plt.subplots(figCount, len(scoreTypes), sharey="row")
    fig.set_size_inches(len(scoreTypes)*5,figCount*4)

    plotFname = "out/compEffScores-valid.pdf"

    plotDataset("xu2015TrainHl60", axArr[0], "Wang 2014\nhuman HL-60\nAs used by Xu 2015", diam=2, addColLabels=True, yLabel="-log2 sgRNA fold change")
    #plotDataset("xu2015TrainKbm7", axArr[1], "Wang 2014\nhuman KBM7\nAs used by Xu 2015", diam=2, addColLabels=True, yLabel="-log2 sgRNA fold change")
    plotDataset("doench2014-Hs", axArr[1], "Doench 2014\nhuman MOLM13, NB4, TF1", diam=2, yLabel="rank-percent", ylim=(0,1.0))
    #plotDataset("doench2014-Mm", axArr[2], "Doench 2014\nmouse EL4", diam=2, yLabel="rank-percent", ylim=(0,1.0))
    plotDataset("chari2015Train293T", axArr[2], "Chari 2015\nhuman 293T", diam=2, yLabel="Mutation Rate", ylim=(0,2.0))
    #plotDataset("chari2015TrainK562", axArr[4], "Chari 2015\nhuman K562", diam=2, yLabel="Mutation Rate")

    plotDataset("xu2015FOX-AR", axArr[3], "Xu 2015 validation\nLNCaP-abl cells, FOX/AR locus\nLentivirus, Western Blot", addColLabels=True)
    plotDataset("xu2015AAVS1",  axArr[4], "Xu 2015 validation\nLNCaP-abl cells, AAVS1 locus\nLentivirus, T7", addColLabels=True)
    axStart = 5
    #chariCells = ["293T", "A549", "HepG2", "K562", "PGP1iPS", "SKNAS", "U2OS"]
    chariCells = [("293T", "Transfection (Lipofect.)"), ("HepG2", "Transfection (Lonza 4D.)")]
    for i in range(0, len(chariCells)):
        cell, transfDesc = chariCells[i]
        dataset = "chari2015Valid_"+cell
        datasetDescs[dataset] = "Chari 2015 %s Validation" % cell
        plotDataset(dataset, axArr[axStart+i], "Chari 2015\nhuman, %s\n%s, sequencing"%(cell, transfDesc))

    plotDataset("doench2014-CD33Exon2", axArr[7], "Doench\nNB4 cells, CDS33 Exon2\nLentivirus", yLabel="sgRNA fold enrichment")
    #plotDataset("doench2014-CD33Exon3", axArr[3], "hg19", "Doench\nNB4 cells\nCD33 Exon3", yLabel="sgRNA fold enrichment")
    #plotDataset("doench2014-CD13Exon10", axArr[4], "hg19", "Doench\nNB4 cells\nCD13 Exon10", yLabel="sgRNA fold enrichment")
    plotDataset("morenoMateos2015", axArr[8], "Moreno-Mateos 2015\nZebrafish RNA injection", diam=3)
    plotDataset("varshney2015", axArr[9], "Varshney 2015\nZebrafish RNA injection")
    plotDataset("gagnon2014", axArr[10], "Gagnon 2014\nZebrafish RNA injection")
    plotDataset("ren2015", axArr[11],  "Ren 2015\nDrosophila injection")
    plotDataset("housden2015", axArr[12],  "Housden 2015\nDrosophila S2R+ cells\nLuciferase-assay")
    plotDataset("farboud2015", axArr[13],  "Farboud 2015\nC. elegans injection")

    #plotDataset("museumT7", axArr[9],  "Concordet\ncell type?, PPP1R12C locus\nelectrop., T7")
    #plotDataset("museumIC50", axArr[10],  "Concordet\ncells?\nelectrop., IC50 assay(name?)")
    plotDataset("schoenig", axArr[14],  "Schoenig\nK562\nLipofection (K2), bGal assay\nbGal: Wefers, PNAS 2013", yLabel="relative rank: 3 (best), 2 or 1", yTicks=[1,2,3])
    #plotDataset("eschstruth", axArr[12],  "Eschstruth\nZebrafish\nInjection", yLabel="relative rank: 3 (best), 2 or 1", yTicks=[1,2,3])
    #plotDataset("concordet2-Hs", axArr[9],  "")
    #plotDataset("concordet2-Mm", axArr[10], "mm9", "")
    #plotDataset("concordet2-Rn", axArr[11], "rn5", "")
    #plotDataset("concordet2", axArr[13],  "Concordet\nhuman/mouse/rat, cellType?\nDelivery?", yLabel="relative rank: 3(best), 2 or 1")

    global datasetDescs

    fig.tight_layout()
    fig.subplots_adjust(left=0.15, top=0.95)
    fig.savefig(plotFname, format = 'pdf')
    fig.savefig(plotFname.replace(".pdf", ".png"))
    print "wrote plot to %s, added .png" % plotFname
    plt.close()

    if scoreCorrFh is not None:
        print "wrote score R summary to %s" % scoreCorrFh.name

def main():
    #"XuData/modFreq.tab" 
    #extendTabAddContext("temp.tab")

    plotLargeScale("out/effScoreComp.tsv")
    plotSmallScale()

main()
