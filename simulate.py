# draw 2,3,4,5... guides with at least eff score 0.1, 0.2, 0.3, 0.4, 0.5 ... from a dataset 
# repeat 100 times
# and show how often we got at least one with a modification frequency > 70%

from annotateOffs import *

def filterGuides(guideData, scoreType, minScore):
    " remove all guides with a given score < minScore "
    filtGuides = []
    for guide in guideData:
        guideName, modFreq, scores = guide
        if scores[scoreType] > minScore:
            filtGuides.append(guide)
    return filtGuides

def runSimulation(fname, reqMinFreq):
    print "Using %s" % fname
    guideData = []
    for row in iterTsvRows(fname):
        scores = dict()
        scores["doench"] = float(row.doench)
        scores["svm"] = float(row.svm)
        scores["ssc"] = float(row.ssc)
        guideData.append( (row.guide, float(row.modFreq), scores) )

    print "min req mod freq: %f" % reqMinFreq
    print

    for sampleSize in [1, 2, 3]:
        print "sample size: %d" % sampleSize
        
        for scoreType in ["svm", "doench", "ssc"]:
            print "score: ", scoreType
            for minScore in range(0, 100, 20):
                minScoreFloat = minScore/100.0
                filtGuides = filterGuides(guideData, scoreType, minScoreFloat)
                #print filtGuides
                okCount = 0
                if len(filtGuides)<10:
                    #print "minScore %f: not enough guides" % minScore
                    continue
                for i in range(0, 1000):
                    # take two guides and check if they have at least one with the req min freq
                    selGuides = random.sample(filtGuides, sampleSize)
                    #print selGuides
                    for g in selGuides:
                        if g[1] > reqMinFreq:
                            okCount +=1
                            break
                succChance = float(okCount)/1000.0
                print "sampleSize: %d, score > %f (%d guides): %0.2f chance" % \
                    (sampleSize, minScoreFloat, len(filtGuides), succChance)
        print 


def main():
    #runSimulation("out/xu2015-compDoenchSsc.tsv", 80)
    runSimulation("out/varshney2015-compDoenchSsc.tsv", 20)

main()
