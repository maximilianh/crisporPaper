# draw 2,3,4,5... guides with at least eff score 0.1, 0.2, 0.3, 0.4, 0.5 ... from a dataset 
# repeat 100 times
# and show how often we got at least one with a modification frequency > 70%

from annotateOffs import *
import numpy as np

def filterGuides(guideData, scoreType, minScore):
    " remove all guides with a given score < minScore "
    filtGuides = []
    for guide in guideData:
        modFreq, scores = guide
        if scores[scoreType] > minScore:
            filtGuides.append(guide)
    return filtGuides

def getScoreRanges(inDir, percentile):
    " for each score, get its range, return as dict scoreName -> (min, max) "
    scoreTypes = ["doench", "ssc", "wangOrig", "chariRaw", "finalGc6"]

    ranges = {}
    for st in scoreTypes:
        ranges[st] = [0.0,0.0]

    # calc all seq scores
    # split them into fname -> list of (koEff, scores)
    fileGuideData = {}
    fileMinReqFreqs = {}
    for fname in glob.glob(inDir+"/*.scores.tab"):
        guideData = []
        koEffs = []
        dataName = basename(fname).split(".")[0]

        for row in iterTsvRows(fname):
            freq = float(row.modFreq)
            rowDict = row._asdict()

            koEff = float(row.modFreq)
            koEffs.append(koEff)

            scores = {}
            for scoreName in row._fields[7:]:
                scores[scoreName] = rowDict[scoreName]
            guideData.append( (koEff, scores) )

        koEffs.sort()
        fileMinReqFreqs[dataName] = koEffs[int(len(koEffs)*percentile)]
        fileGuideData[dataName] = guideData
    return fileGuideData, fileMinReqFreqs

    #fileGuideData = {}
    #fileMinReqFreqs = {}
    #for fname in glob.glob(inDir+"/*.ext.tab"):
        #guideData = []
        #koEffs = []
        #for row in iterTsvRows(fname):
            #seq = row.extSeq
            #scores = seqScores[seq]
            #koEff = seqKoEffs[seq]
            #guideData.append( (seq, koEff, scores) )
            #koEffs.append(koEff)

        #dataName = basename(fname).split(".")[0]
        #koEffs.sort(reverse=True)
        #fileMinReqFreqs[dataName] = koEffs[int(len(koEffs)*percentile)]
        #fileGuideData[dataName] = guideData
    #return ranges, guideData
    #return fileGuideData, fileMinReqFreqs

def runSimulation(inDir, percentile):
    #print "Using %s" % fname
    #guideData = []
    #for row in iterTsvRows(fname):
    #    scores = dict()
    #    scores["doench"] = float(row.doench)
    #    scores["svm"] = float(row.svm)
    #    scores["ssc"] = float(row.ssc)
    #    guideData.append( (row.guide, float(row.modFreq), scores) )

    guideDatas, reqMinFreqs  = getScoreRanges(inDir, percentile)

    print

    # print realRanges
    scoreRanges = {'doench': [0.0, 1.0], 'wangOrig': [0.0, 1.0], 'chariRaw': [-3.0, +3.0], 'ssc': [-1.5, 1.5], 'finalGc6': [0.0, 6.0], 'fusi' : (0.0,1.0)}

    for dataset, guideData in guideDatas.iteritems():
        reqMinFreq = reqMinFreqs[dataset]
        if not "xu2015Train" in dataset:
            continue
        print
        print "dataset %s, min req mod freq: %f" % (dataset, reqMinFreq)

        for sampleSize in [1, 2, 3]:
            print "sample size: %d" % sampleSize
            
            for scoreType in ["wangOrig", "doench", "ssc", "chariRaw", "fusi", "finalGc6"]:
                print "score: ", scoreType
                minScore, maxScore = scoreRanges[scoreType]
                for minScoreFloat in np.linspace(minScore, maxScore, 20):
                    filtGuides = filterGuides(guideData, scoreType, minScoreFloat)
                    #print filtGuides
                    okCount = 0
                    okSum = 0
                    if len(filtGuides)<10:
                        #print "minScore %f: not enough guides" % minScore
                        continue
                    for i in range(0, 1000):
                        # take two (or more) guides and check if they have at least one with the req min freq
                        selGuides = random.sample(filtGuides, sampleSize)
                        #print selGuides
                        foundOk = 0
                        for g in selGuides:
                            if g[0] > reqMinFreq:
                                foundOk += 1
                        if foundOk > 0:
                            okCount += 1
                        okSum += foundOk
                    avgOk = float(okSum) / 1000
                    succChance = float(okCount)/1000.0
                    print "score > %f (%d guides): %0.2f chance, avg good guides %0.2f" % \
                        (minScoreFloat, len(filtGuides), succChance, avgOk)
            print 


def main():
    #runSimulation("out/xu2015-compDoenchSsc.tsv", 80)
    #runSimulation("out/varshney2015-compDoenchSsc.tsv", 20)
    # out/xu2015Train-compEffData.tsv
    runSimulation("effData", 0.8)

main()
