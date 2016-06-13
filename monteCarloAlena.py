# determine how well some score predicts performance in the alena dataset
# on the level of the locus, not overall
from annotateOffs import *
from collections import defaultdict
import random
import operator


def parseGuides(fname):
    " parse scores.tab, return dict locusName -> (score1, score2, modFreq) "
    guideCount = 0
    scoreData = defaultdict(list)
    for row in iterTsvRows(fname):
        setName = row.guide.split("_")[0]
        predScore1 = float(row.crisprScan)
        #predScore1 = float(row.doench)
        #predScore2 = float(row.doench)
        predScore2 = float(row.fusi)
        scoreData[setName].append((predScore1, predScore2, float(row.modFreq)))
        guideCount += 1
    return scoreData, guideCount

def getTopFreqs(scoreData):
    """
    return dict with locusNames -> list of top frequencies
    """
    locusTopFreqs = dict()
    for locusName, locusGuides in scoreData.items():
        modFreqs = [modFreq for predScore1,predScore2,modFreq in locusGuides]
        modFreqs.sort()
        topFreqs = modFreqs[-2:]
        locusTopFreqs[locusName] = topFreqs
    return locusTopFreqs

def countSuccesses(scoreData, locusTopFreqs, predFunc):
    """ for each locus, identify best guides using selection function and
    count how often we find indeed a top guide
    topScores is a dict locus -> set of top scores
    """
    successCount = 0
    for locusName, locusGuides in scoreData.items():
        # sort by pred score and get data of top three pred scores
        # get the best two mod freqs
        topFreqs= locusTopFreqs[locusName]
        topGuides = predFunc(locusGuides)
        # check if any of the top guides has a top mod freq
        success = False
        for topGuide in topGuides:
            predScore1, predScore2,modFreq = topGuide
            if modFreq in topFreqs:
                success = True
                break
        if success:
            successCount += 1
    return successCount
        #print "Locus and guideCount", locusName, len(locusGuides), "best pred:", topGuides, "highest measured mod. freqs:", topFreqs, "prediction worked?", success
    #print "%d guides, "%guideCount, "number of loci", len(scoreData), "success count", successCount


def getBestCombGuides(guides):
    " return a set of guides predicted by score1 and score2 "
    guides.sort()
    topGuides = []
    # take top2 based on primary score
    topGuides.extend(guides[-1:])
    guides.sort(key=operator.itemgetter(1))
    topGuides.extend(guides[-1:])
    assert(len(topGuides)==2)
    return topGuides

def getBestTwoPredGuides(guides):
    " return a set of guides predicted by score1 "
    guides.sort()
    topGuides = []
    # take top2 based on primary score
    topGuides.extend(guides[-2:])
    return topGuides

def getBestPredGuide(guides):
    " return a set of guides predicted by score1 "
    guides.sort()
    topGuides = []
    # take top2 based on primary score
    topGuides.extend(guides[-1:])
    return topGuides

def getRndGuides(guides):
    " return a set of four random guides "
    random.shuffle(guides)
    return guides[:2]

print "Shkumatava lab zebrafish dataset"
scoreData, guideCount = parseGuides("effData/alenaAll.scores.tab")
topFreqs = getTopFreqs(scoreData)
print "Strategy: pick 2 guides with highest Moreno-Mateos score"
successCount = countSuccesses(scoreData, topFreqs, getBestTwoPredGuides)
#successCount = countSuccesses(scoreData, topFreqs, getBestCombGuides)
print "success: ", successCount, "out of", len(scoreData), "loci"

n = 100000
rndSuccCount = 0
for i in range(0, n):
    rndCount = countSuccesses(scoreData, topFreqs, getRndGuides)
    if rndCount >= successCount:
        rndSuccCount +=1
print "pVal", rndSuccCount / float(n)

# copied for the Teboul dataset
print "Teboul in-vivo dataset"
scoreData, guideCount = parseGuides("effData/teboulVivo_mm9.scores.tab")
topFreqs = getTopFreqs(scoreData)
print "Strategy: pick 1 guide with highest Moreno-Mateos score"
successCount = countSuccesses(scoreData, topFreqs, getBestPredGuide)
#successCount = countSuccesses(scoreData, topFreqs, getBestCombGuides)
print "success: ", successCount, "out of", len(scoreData), "loci"

n = 100000
rndSuccCount = 0
for i in range(0, n):
    rndCount = countSuccesses(scoreData, topFreqs, getRndGuides)
    if rndCount >= successCount:
        rndSuccCount +=1
print "pVal", rndSuccCount / float(n)

