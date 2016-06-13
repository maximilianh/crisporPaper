# plot the binary classification metrics

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
import matplotlib.pyplot as plt
import numpy as np
from annotateOffs import *
from collections import defaultdict
import random

scoreNames = ['doench', 'ssc', 'crisprScan', 'wangOrig', 'chariRank', "wuCrispr", 'fusi', 'finalGc6', 'finalGg']
# doench = regression
# ssc = regression
# crisprScan = regression
# wang = SVM
# chari = SVM

dataSubs = {
     'varshney2015': ("Zebrafish", "Injection", "Sequencing"),
     'ren2015': ("Drosophila", "Injection", "Sequencing"),
     'xu2015TrainHl60': ("KBM7/HL60", "Lentivir.", "KO"),
     'gagnon2014': ("Zebrafish", "Injection", "Sequencing"),
     'chari2015Train': ("293T", "Lentivir.", "Lib-on-Lib KO"),
     'chari2015Valid_293T':("293T", "Transfection", "Sequencing"),
     'doench2014-Hs': ("MOLM13/NB4/TF1", "Lentivir.",  'KO'),
     #'museumIC50':  ("?", "?", "?"),
     'xu2015AAVS1': ("LNCap-abl", "Lentivir.", "West.Blot"),
     'xu2015FOX-AR': ("LNCap-abl", "Lentivir.", "T7 Endonucl."),
     'schoenig': ("K562", "betaGal-assay", "betaGal"),
     'farboud2015' : ("Zebrafish", "Injection", "Sequencing"),
     #'eschstruth' : ("Zebrafish", "Injection", "T7 Endonucl."),
     'morenoMateos2015' : ("Zebrafish", "Injection", "Sequencing"),
     'alenaAll' : ("Zebrafish", "Injection", "Sanger Seq"),
     'hart2016-Hct1162lib1Avg' : ("Hct116 2", "Lentivir.", "Sequencing"),
     'ghandi2016_ci2' : ("Ciona", "Electroporation", "Sequencing"),
     'teboulVivo_mm9' : ("Mouse", "Injection", "Mutant Embryos"),
     'concordet2' : ("U2OS/MEF/C6", "Electrop.", "T7 Endon."),
     #'housden2015' : ("Dros. S2R+", "Transfection", "Lucif.")
    }
topDatasets = [
     'xu2015TrainHl60',
     'doench2014-Hs',
     'chari2015Train',
     'farboud2015',
     'ren2015',
     #'housden2015',
    'hart2016-Hct1162lib1Avg',
    'ghandi2016_ci2',
    ]

middleDatasets = [
    #'xu2015AAVS1',
    #'xu2015FOX-AR',
    #'chari2015Valid_293T',
     'morenoMateos2015',
    "varshney2015",
    "gagnon2014"
]

#buttomDatasets = [
     #'ren2015',
     #'farboud2015',
     #'gagnon2014',
     #'varshney2015',
#]

scoreDescs = {
    "wangOrig" : "Wang",
    "doench" : "Doench",
    "ssc" : "Xu (Wang)",
    "chariRank" : "Chari Rank",
    "crisprScan" : "Moreno-Mateos",
    "fusi" : "Fusi/Doench",
    "chariRaw" : "Chari",
    "finalGc6" : "Ren: 3'GC>4",
    #"drsc" : "Housden",
    "wuCrispr" : "Wong",
    #"finalGc2" : "Farboud-like, last 2 bp GC",
    "finalGg" : "Farboud: -GG",
}
def parseData(fname):
    """ return dict of scoreType -> tuple of (recallList, precisionList, f1List) 
    (one element per dataset) and the names of the datasets"""
    scoreDict = defaultdict(dict)
    dataMax = {}
    dataCountInfo = {}
    for row in iterTsvRows(fname):
        if row.classifierName.startswith("DecTree"):
            continue
        rec = int(float(row.recall)*100)
        prec = int(float(row.precision)*100)
        f1 = int(float(row.f1)*100)
        scoreDict[row.classifierName][row.dataset] = (rec, prec, f1)
        dataMax[row.dataset] = max(dataMax.get(row.dataset, 0), f1)
        dataCountInfo[row.dataset] = (int(float(row.size)), int(float(row.posCount)))

    scoreDict = dict(scoreDict)
    #print scoreDict

    # sort dataNames by f1 value
    dataMaxes = dataMax.items()
    dataMaxes.sort(key=operator.itemgetter(-1), reverse=True)
    dataNames = [x for x,y in dataMaxes]

    # but put some first
    newDataNames = []
    newDataNames.extend(topDatasets)
    newDataNames.extend(middleDatasets)
    for dn in dataNames:
        if dn not in topDatasets and dn not in middleDatasets:
            newDataNames.append(dn)
        #if not dn in topDatasets and not dn in middleDatasets:
            #newDataNames.append(dn)
    #for dn in dataNames:
        #if dn in middleDatasets:
            #newDataNames.append(dn)
    #for dn in dataNames:
        #if dn in topDatasets:
            #newDataNames.append(dn)
    dataNames = list(reversed(newDataNames))

    # transform to dict className -> list of scores
    scores = dict()
    for scoreName in scoreNames:
        dataVals = scoreDict[scoreName]
        recList = []
        precList = []
        f1List = []
        for dataName in dataNames:
            rec, prec, f1 = dataVals[dataName]
            recList.append(rec)
            precList.append(prec)
            f1List.append(f1)
        scores[scoreName] = (recList, precList, f1List)

    return scores, dataNames, dataCountInfo

def plot(scores, dataNames, dataCountInfo, outfname):
    " "
    plt.figure(figsize=(8,12))
    plt.rcParams['ytick.major.pad']='8'

    plots = []
    colors = list(reversed(["blue", "red", "orange", "magenta", "orange", "red", "blue", "black", "black"]))
    markers = list(reversed(["^",   "^",   "^",      "o",       "o",      "s",    "s",      "x",      "+" ]))

    fig, axArr = plt.subplots(1, 2, sharey=True)

    for scoreName in scoreNames:
        dataTuple = scores[scoreName]
        recVals, precVals, f1Vals = dataTuple
        yPosList = range(0, len(f1Vals))
        yPosList = [y-random.random()*0.25 for y in yPosList]
        col = colors.pop()
        marker = markers.pop()
        alpha = 0.7
        plot = axArr[0].scatter(precVals, yPosList, alpha=alpha, s=30, color=col, marker=marker)
        plot = axArr[1].scatter(recVals, yPosList, alpha=alpha, s=30, color=col, marker=marker)
        plots.append(plot)

    #plots = []
    axArr[1].legend(plots,
           [scoreDescs[x] for x in scoreNames],
           labelspacing=0,
           #bbox_to_anchor = (0,0,1,1),
           #bbox_transform = plt.gcf().transFigure ,
           scatterpoints=1,
           loc='upper right',
           #ncol=len(scoreNames),
           fontsize=8)

    plt.setp(axArr[1].get_yticklabels(), visible=False)
    axArr[0].set_ylim(-1,len(dataNames))
    axArr[0].set_yticks(range(0, len(dataNames)))
    axArr[0].set_yticklabels([mainDataDescs[x] for x in dataNames])
    axArr[0].set_xlim(-5,105)
    axArr[1].set_xlim(-5,105)
    ls = ":"
    lw = 0.4
    for i in [0,1]:
        axArr[i].axhline(len(dataNames)-len(topDatasets)-0.7, ls=ls, color="k", lw=lw)
        axArr[i].axhline(len(dataNames)-len(topDatasets)-len(middleDatasets)-0.7, ls=ls, color="k", lw=lw)
    axArr[0].set_xlabel("Precision")
    axArr[1].set_xlabel("Recall")
    [tick.label.set_fontsize(10) for tick in axArr[0].yaxis.get_major_ticks()]

    fig.tight_layout()
    #for tl in  axArr[0].get_yticklabels():
        #print tl.get_position()
    for y in range(0, len(dataNames)):
        dataName = dataNames[y]
        dataSubStr = " - ".join(dataSubs[dataName])
        axArr[0].annotate(dataSubStr, xy=(0,0), ha="right", size="8", xytext=(-8, y-0.4))
        size, posCount = dataCountInfo[dataName]
        axArr[0].annotate("%d guides, %d positives" % (size, posCount), xy=(0,0), ha="right", size="8", xytext=(-8, y-0.68))

    fig.subplots_adjust(wspace=0.0001)
    #fig = plt.gcf()
    fig.subplots_adjust(left=0.27)
    plt.savefig(outfname)
    print "wrote %s" % outfname

    outfname = outfname.replace("pdf", "png")
    plt.savefig(outfname)
    print "wrote %s" % outfname
    plt.close()

def main():
    scores, dataNames, dataCountInfo = parseData("out/binClassMetrics.tsv")
    plot(scores, dataNames, dataCountInfo, "out/precRecall.pdf")


main()
