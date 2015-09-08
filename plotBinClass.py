# plot the binary classification metrics

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
import matplotlib.pyplot as plt
import numpy as np
from annotateOffs import *
from collections import defaultdict

scoreNames = ['finalGg', 'finalGc6', 'svm', 'doench', 'ssc', 'chariRank']

dataDescs = {
     'varshney2015': "Varshney Zebrafish",
     'ren2015': "Ren Drosophila",
     'xu2015Train': "Wang/Xu KO",
     'gagnon2014': "Gagnon Zebrafish",
     'chari2015Train':"Chari Training",
     'chari2015Valid_293T':"Chari Validation",
     'doench2014-Hs': 'Doench Training',
     'museumIC50':  "Concordet IC50",
     'xu2015AAVS1': "Xu Validation AAVS1",
     'xu2015FOX-AR': "Xu Validation FOX/AR",
     'schoenig': "Schoenig LacZ",
     'farboud2015' : "Farboud C.elegans",
     'eschstruth' : "Eschstruth Zebrafish"
     }

dataSubs = {
     'varshney2015': ("Eggs", "Injection", "Sequencing"),
     'ren2015': ("Eggs", "Injection", "Sequencing"),
     'xu2015Train': ("KBM7/HL60", "Lentivir.", "KO"),
     'gagnon2014': ("Eggs", "Injection", "Sequencing"),
     'chari2015Train': ("293T", "Lentivir.", "Lib-on-Lib KO"),
     'chari2015Valid_293T':("293T", "Transfection", "?"),
     'doench2014-Hs': ("MOLM13/NB4/TF1", "Lentivir.",  'KO'),
     'museumIC50':  ("?", "?", "?"),
     'xu2015AAVS1': ("LNCap-abl", "Lentivir.", "West.Blot"),
     'xu2015FOX-AR': ("LNCap-abl", "Lentivir.", "T7"),
     'schoenig': ("K562", "betaGal-assay", "?"),
     'farboud2015' : ("Eggs", "Injection", "Sequencing"),
     'eschstruth' : ("Eggs", "Injection", "T7")
    }
topDatasets = [
     'xu2015Train',
     'chari2015Train',
     'doench2014-Hs'
    ]

middleDatasets = [
     'ren2015',
     'farboud2015',
     'gagnon2014',
     'varshney2015',
]

scoreDescs = {
    "svm" : "Wang Score",
    "doench" : "Doench Score",
    "ssc" : "Xu Score",
    "chariRank" : "Chari Score",
    "finalGc6" : "Ren: 3' 6bp GC>4",
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

    # sort dataNames by f1 value
    dataMaxes = dataMax.items()
    dataMaxes.sort(key=operator.itemgetter(1), reverse=False)
    dataNames = [x for x,y in dataMaxes]

    # and put some last
    newDataNames = []
    for dn in dataNames:
        if not dn in topDatasets and not dn in middleDatasets:
            newDataNames.append(dn)
    for dn in dataNames:
        if dn in middleDatasets:
            newDataNames.append(dn)
    for dn in dataNames:
        if dn in topDatasets:
            newDataNames.append(dn)
    dataNames = newDataNames

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
    plt.rcParams['ytick.major.pad']='8'

    plots = []
    colors = ["black", "violet", "red", "green", "orange", "black", "blue"]
    markers = list(reversed(["x", "+", "o", "s", ">", "<", "^"]))

    fig, axArr = plt.subplots(1, 2, sharey=True)

    for scoreName in scoreNames:
        dataTuple = scores[scoreName]
        recVals, precVals, f1Vals = dataTuple
        xPosList = range(0, len(f1Vals))
        col = colors.pop()
        marker = markers.pop()
        plot = axArr[0].scatter(recVals, xPosList, alpha=0.6, s=30, color=col, marker=marker)
        plot = axArr[1].scatter(precVals, xPosList, alpha=0.6, s=30, color=col, marker=marker)
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
           fontsize=9)

    plt.setp(axArr[1].get_yticklabels(), visible=False)
    axArr[0].set_ylim(-1,len(dataNames))
    axArr[0].set_yticks(range(0, len(dataNames)))
    axArr[0].set_yticklabels([dataDescs[x] for x in dataNames])
    axArr[0].set_xlim(-5,105)
    axArr[1].set_xlim(-5,105)
    ls = ":"
    axArr[0].axhline(len(dataNames)-len(topDatasets)-0.7, ls=ls, color="k", lw=0.1)
    axArr[1].axhline(len(dataNames)-len(topDatasets)-0.7, ls=ls, color="k", lw=0.1)
    axArr[0].axhline(len(dataNames)-len(topDatasets)-len(middleDatasets)-0.7, ls=ls, color="k", lw=0.1)
    axArr[1].axhline(len(dataNames)-len(topDatasets)-len(middleDatasets)-0.7, ls=ls, color="k", lw=0.1)
    #axArr[0].xlabel("F1 Score (harmonic mean of precision and recall)")
    axArr[0].set_xlabel("Recall")
    axArr[1].set_xlabel("Precision")
    [tick.label.set_fontsize(10) for tick in axArr[0].yaxis.get_major_ticks()]

    fig.tight_layout()
    for tl in  axArr[0].get_yticklabels():
        print tl.get_position()
    for y in range(0, len(dataNames)):
        dataName = dataNames[y]
        dataSubStr = " - ".join(dataSubs[dataName])
        axArr[0].annotate(dataSubStr, xy=(0,0), ha="right", size="8", xytext=(-8, y-0.37))
        size, posCount = dataCountInfo[dataName]
        axArr[0].annotate("%d guides, %d positives" % (size, posCount), xy=(0,0), ha="right", size="8", xytext=(-8, y-0.62))

    fig.subplots_adjust(wspace=0.0001)
    plt.savefig(outfname)
    print "wrote %s" % outfname

    outfname = outfname.replace("pdf", "png")
    plt.savefig(outfname)
    print "wrote %s" % outfname
    plt.close()

def main():
    scores, dataNames, dataCountInfo = parseData("out/binClassMetrics.tsv")
    plot(scores, dataNames, dataCountInfo, "out/binClass.pdf")


main()
