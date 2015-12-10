# for each dataset, get the top25% and the bottom25% and assign labels 1.0 and 0.0
# then get all scores for these sequences and write scores and labels to a file precRec/dataset-scoreType.tab
# then use the R ROCR package and write the x and y values to a file
# precRec/dataset-scoreType.precRec.tab 
# finally plot all these curves with matplotlib

import os, logging, glob
from annotateOffs import *
from collections import defaultdict
from os.path import *

import matplotlib
#matplotlib.use('Agg')
#matplotlib.rcParams['pdf.fonttype'] = 42
import matplotlib.pyplot as plt

selDatasets = ["doench2014-Hs", "morenoMateos2015", "chari2015Train", "xu2015TrainHl60", "chariEval"]

scoreDescs = {
    "wang" : "Wang Score2",
    "wangOrig" : "Wang Score",
    "doench" : "Doench Score",
    "ssc" : "Xu Score",
    "chariRaw" : "Chari Score",
    "chariRank" : "Chari Rank Score",
    "chariRaw" : "Chari Score",
    "finalGc6" : "Ren: last 6 bp GC>4",
    #"finalGc2" : "Farboud-like, last 2 bp GC",
    "finalGg" : "Farboud: ends with GG",
    "myScore" : "Max Score",
    "crisprScan" : "Moreno-Matos Score",
    "fusi" : "Fusi/Doench Score",
    "drsc" : "Housden Score",
    "wuCrispr" : "Wong Score",
    "oof"      : "Bae Out-of-Frame Score"
}

# doench = regression
def getQuartiles(freqs):
    """ given a seq -> (name, freq) dict return the sequences with top/bottom 25% freq """
    freqNames = []
    for seq, nameFreq in freqs.iteritems():
        name, freq = nameFreq
        freqNames.append( (freq, seq) )
    freqNames.sort()

    topQuartStart = int(len(freqNames)*0.75)
    bottomQuartEnd = int(len(freqNames)*0.25)
    topQuart = freqNames[topQuartStart:]
    bottomQuart = freqNames[:bottomQuartEnd]

    topSeqs = [y for x,y in topQuart]
    bottomSeqs = [y for x,y in bottomQuart]
    return topSeqs, bottomSeqs

def getScores(seqScores, seqs):
    """ given a dict seq -> seqType -> score and a list of seqs, return a dict 
    scoreType -> list of score (float, same order as in seqs)
    """
    nameToScores = {}
    for scoreName in allScoreNames:
        scores = []
        for s in seqs:
            scores.append(seqScores[s][scoreName])
        nameToScores[scoreName] = scores

    return nameToScores

def readData(datasetName):
    " return the 1/0 labels and a dict with scoreName -> list of scores "
    seqScores, freqs = parseEffScores(datasetName)
    topSeqs, bottomSeqs = getQuartiles(freqs)

    labels = [1]*len(topSeqs)
    labels.extend( [0]*len(bottomSeqs) )
    allSeqs = topSeqs
    allSeqs.extend(bottomSeqs)
    allScores = getScores(seqScores, allSeqs)
    #assert(len(labels)==len(allScores.values()[0]))

    return labels, allScores, allSeqs

def writeSeqLabelsScores(allSeqs, labels, allScores, outDir, datasetName):
    for scoreName, scores in allScores.iteritems():
        assert(len(allSeqs)==len(scores)==len(labels))
        outFname = join(outDir, "%s_%s.tab" % (datasetName, scoreName))
        ofh = open(outFname, "w")
        ofh.write("seq\tlabel\tscore\n")
        for seq, label, score in zip(allSeqs, labels, scores):
            ofh.write("%s\t%d\t%f\n" % (seq, label, score))
        ofh.close()
        print "wrote %s" % ofh.name

def writeRTables(outDir, datasetName):
    " write R dataframes to outDir/datasetName-scoreName.tab "
    allFound = True
    for scoreName in allScoreNames:
        outFname = join(outDir, "%s_%s.tab" % (datasetName, scoreName))
        if not isfile(outFname):
            allFound = False
            
    if allFound:
        print "all scoring files already exist for dataset %s" % datasetName
        return

    labels, allScores, allSeqs = readData(datasetName)
    writeSeqLabelsScores(allSeqs, labels, allScores, outDir, datasetName)

def parseCoords(inDir, datasetName):
    """ parse all files in inDir/datasetName*.coords.tab and return as dict
        scoreName -> (xList, yList)
    """
    ret = {}
    fnames = glob.glob(join(inDir, "%s_*.coords.tab" % datasetName))
    if len(fnames)==0:
        raise Exception("No filenames for dataset %s" % datasetName)

    for fname in fnames:
        xList, yList = [], []
        for line in open(fname):
            x, y = line.strip().split()
            if y=="NA":
                continue
            xList.append(float(x))
            yList.append(float(y))
        
        scoreName = basename(fname).split('.')[0].split('_')[-1]
        ret[scoreName] = (xList, yList)
    return ret

def runR(dataDir):
    " run R on all files in dataDir and generate the x/y coordinates "
    for inFname in glob.glob(join(dataDir, "*.tab")):
        if "coords.tab" in inFname:
            continue
        outFname = inFname.replace(".tab", ".coords.tab")
        if isfile(outFname):
            continue
        cmd = "Rscript precRecCurve.R %s %s"% (inFname, outFname)
        assert(os.system(cmd)==0)

def plotCoords(coords, axis):
    " plot the x-y coordinates onto the axis "
    plots = []
    lineStyles = ['-', '--', ':', "-", "--", ":", "-", "--", "-", "-"]
    selColors = ["blue", "red", "orange", "magenta", "orange", "black", "orange", "black", "black", "black"]

    legendNames = []
    for i, scoreName in enumerate(allScoreNames):
        if scoreName.startswith("final"):
            continue
        if scoreName=="drsc":
            continue
        legendNames.append(scoreDescs[scoreName])
        xyLists = coords[scoreName]
        xList, yList = xyLists
        color = selColors[i]
        lineStyle = lineStyles[i]
        xList, yList = xyLists
        plot, = axis.plot(xList, yList, color=color, ls=lineStyle)
        plots.append(plot)
    axis.set_xlim(0,1.0)
    axis.set_ylim(0,1.0)
    #print scoreName, xList, yList
    #plot = axArr[0].scatter(xList, yList)
    #scoreNames.append(scoreName)
    return plots, legendNames

def makeChariEval():
    " create a special dataset chariEval based on S1 of Chari et al "
    #labels, allScores, allSeqs = readData('chari2015TrainK562')
    seqScores, freqs = parseEffScores('chari2015Train293T')
    topSeqs = open("effData/chari2015Train/chari-S1-high-CasSp.txt").read().splitlines()
    bottomSeqs = open("effData/chari2015Train/chari-S1-low-CasSp.txt").read().splitlines()

    labels = [1]*len(topSeqs)
    labels.extend([0]*len(bottomSeqs))

    evalSeqs = list(tuple(topSeqs)) # copy of list
    evalSeqs.extend(bottomSeqs)

    #wuCrisprFilter = False
    wuCrisprFilter = True
    if wuCrisprFilter:
        # keep only sequences where the wuCrispr score is not 0
        filtSeqs = []
        for seq in evalSeqs:
            scores = seqScores[seq]
            if scores["wuCrispr"]==0.0:
                print "skip", seq
                continue
            filtSeqs.append(seq)

        print "XXX", filtSeqs
        # make labels for them
        filtLabels = []
        for seq in filtSeqs:
            if seq in topSeqs:
                #print "Yes", seq
                filtLabels.append(1)
            elif seq in bottomSeqs:
                #print "No", seq
                filtLabels.append(0)
            else:
                assert(False)

        print "Chari eval: after filtering, having %d seqs, %d labels" % (len(filtSeqs), len(filtLabels))
        labels = filtLabels
        evalSeqs = filtSeqs

    # and get all scores for them
    evalScores = defaultdict(list)
    for seq in evalSeqs:
        scores = seqScores[seq]
        for scoreName, score in scores.iteritems():
            evalScores[scoreName].append(score)

    print "Chari eval: got %d seqs, %d labels, %d scores" % (len(evalSeqs), len(labels), len(evalScores.values()[0]))
    print labels
    return labels, evalScores, evalSeqs


def main():
    datasetNames = selDatasets
    dataDir = "precRecData"

    fig, axes = plt.subplots(1, len(datasetNames))
    fig.set_size_inches(5*len(datasetNames), 5)

    if not isdir(dataDir):
        os.mkdir(dataDir)
    for datasetName in datasetNames:
        #print "PROCESSING DATASET", datasetName
        if datasetName=="chariEval":
            labels, allScores, allSeqs = makeChariEval()
            writeSeqLabelsScores(allSeqs, labels, allScores, dataDir, datasetName)
        else:
            writeRTables(dataDir, datasetName)
        runR(dataDir)

    for datasetName, axis in zip(datasetNames, axes):
        coords = parseCoords(dataDir, datasetName)
        plots, scoreNames = plotCoords(coords, axis)
        axis.set_title(datasetDescs[datasetName])
        axis.legend(plots, scoreNames, 'lower right', labelspacing=0, fontsize=10, ncol=2)

    plt.tight_layout()
    outfname = "temp.pdf"
    plt.savefig(outfname)
    print "wrote %s" % outfname
    outfname = outfname.replace("pdf", "png")
    plt.savefig(outfname)
    print "wrote %s" % outfname
    plt.close()


main()




