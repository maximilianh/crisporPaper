# plot the binary classification metrics

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
import matplotlib.pyplot as plt
import numpy as np
from annotateOffs import *
from collections import defaultdict
import random

scoreNames = ['doench', 'ssc', 'crisprScan', 'wangOrig', 'chariRank', 'fusi', "drsc", 'finalGc6', 'finalGg']
# doench = regression
# ssc = regression
# crisprScan = regression
# wang = SVM
# chari = SVM

dataDescs = {
     'varshney2015': "Varshney Zebrafish",
     'ren2015': "Ren Drosophila",
     'xu2015TrainHl60': "Wang/Xu KO Training",
     'gagnon2014': "Gagnon Zebrafish",
     'chari2015Train':"Chari Training",
     'chari2015Valid_293T':"Chari Validation",
     'doench2014-Hs': 'Doench Training',
     'museumIC50':  "Concordet IC50",
     'xu2015AAVS1': "Xu Validation AAVS1",
     'xu2015FOX-AR': "Xu Validation FOX/AR",
     'schoenig': u'Sch\u00F6nig LacZ',
     'farboud2015' : "Farboud C.elegans",
     'eschstruth' : "Eschstruth Zebrafish        ",
     'morenoMateos2015' : "CrisprScan Training",
     'alenaAll' : "Shkumatava Dataset",
     'housden2015' : "Housden Dros. Training",
     }

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
     'xu2015FOX-AR': ("LNCap-abl", "Lentivir.", "T7"),
     'schoenig': ("K562", "betaGal-assay", "betaGal"),
     'farboud2015' : ("Zebrafish", "Injection", "Sequencing"),
     'eschstruth' : ("Zebrafish", "Injection", "T7"),
     'morenoMateos2015' : ("Zebrafish", "Injection", "Sequencing"),
     'alenaAll' : ("Zebrafish", "Injection", "Sanger Seq."),
     'housden2015' : ("Dros. S2R+", "Transfection", "Lucif.")
    }
topDatasets = [
     'xu2015TrainHl60',
     'doench2014-Hs',
     'chari2015Train',
     'farboud2015',
     'ren2015',
     'housden2015',
     'morenoMateos2015'
    ]

middleDatasets = [
    #'xu2015AAVS1',
    #'xu2015FOX-AR',
    #'chari2015Valid_293T',
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
    "crisprScan" : "CrisprScan",
    "fusi" : "Fusi (Doench)",
    "chariRaw" : "Chari",
    "finalGc6" : "Ren: 3'GC>4",
    "drsc" : "Housden",
    #"finalGc2" : "Farboud-like, last 2 bp GC",
    "finalGg" : "Farboud: -GG",
}
def parseData(fname):
    """ return dict of scoreType -> dict dataName of (recallList, precisionList, f1List) 
    (one element per dataset) and the names of the datasets"""
    scoreDict = defaultdict(dict)
    dataMax = {}
    dataCountInfo = {}
    for row in iterTsvRows(fname):
        if row.classifierName.startswith("DecTree"):
            continue
        acc = int(float(row.bestXAcc)*100)
        scoreDict[row.classifierName][row.dataset] = acc
        dataMax[row.dataset] = max(dataMax.get(row.dataset, 0), acc)
        dataCountInfo[row.dataset] = (int(float(row.size)), int(float(row.bestXPredCount)))

    scoreDict = dict(scoreDict)
    #print scoreDict

    # sort dataNames by f1 value
    dataMaxes = dataMax.items()
    dataMaxes.sort(key=operator.itemgetter(-1), reverse=True)
    dataNames = [x for x,y in dataMaxes]

    return scoreDict, dataNames, dataCountInfo

def plot(scores, dataNames, dataCountInfo, outfname):
    " "
    plt.figure(figsize=(5,3))
    plt.rcParams['ytick.major.pad']='8'

    plots = []
    colors = list(reversed(["blue", "red", "orange", "magenta", "orange", "grey", "orange", "black", "black"]))
    markers = list(reversed(["^", "^", "^", "o", "o",  "s", "s", "x", "+"]))

    for scoreName in scoreNames:
        dataDict = scores[scoreName]
        precList = []
        for dataName in dataNames:
            precList.append(dataDict[dataName])

        #recVals, precVals, f1Vals = dataTuple
        yPosList = range(0, len(precList))
        yPosList = [y-random.random()*0.1 for y in yPosList]
        col = colors.pop()
        marker = markers.pop()
        alpha = 0.7
        plot = plt.scatter(precList, yPosList, alpha=alpha, s=30, color=col, marker=marker)
        plots.append(plot)

    lgd = plt.legend(plots,
           [scoreDescs[x] for x in scoreNames],
           labelspacing=0,
           bbox_to_anchor = ((1.43, 1.02)),
           #bbox_transform = plt.gcf().transFigure ,
           scatterpoints=1,
           loc='upper right',
           #ncol=len(scoreNames),
           fontsize=8)
    artists = [lgd]

    #plt.setp(axArr[1].get_yticklabels(), visible=False)
    gca = plt.gca()
    gca.set_ylim(-1,len(dataNames))
    gca.set_yticks(range(0, len(dataNames)))
    gca.set_yticklabels([dataDescs[x] for x in dataNames])
    gca.set_xlim(-5,105)
    gca.set_xlabel("Accuracy, in %")
    ls = ":"
    lw = 0.4
    [tick.label.set_fontsize(10) for tick in gca.yaxis.get_major_ticks()]

    plt.tight_layout()
    for y in range(0, len(dataNames)):
        dataName = dataNames[y]
        dataSubStr = " - ".join(dataSubs[dataName])
        annot = plt.annotate(dataSubStr, xy=(0,0), ha="right", size="8", xytext=(-8, y-0.4))
        artists.append(annot)

        size, posCount = dataCountInfo[dataName]
        annot = plt.annotate("%d guides, %d positives" % (size, posCount), xy=(0,0), ha="right", size="8", xytext=(-8, y-0.68))
        artists.append(annot)

    plt.savefig(outfname, bbox_extra_artists=artists, bbox_inches='tight')
    print "wrote %s" % outfname

    outfname = outfname.replace("pdf", "png")
    plt.savefig(outfname, bbox_extra_artists=artists, bbox_inches='tight')
    print "wrote %s" % outfname
    plt.close()

def main():
    scores, dataNames, dataCountInfo = parseData("out/binClassMetrics.tsv")

    plotDataNames = ["schoenig", "eschstruth", "alenaAll"]
    plot(scores, plotDataNames, dataCountInfo, "out/accuracy.pdf")


main()
