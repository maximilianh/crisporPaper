from annotateOffs import *
from collections import defaultdict
import random

realData = defaultdict(list)
scoreData = defaultdict(list)
for row in iterTsvRows("effData/schoenig.scores.tab"):
    #if row.guide.startswith("cdk"):
        #continue
    #if row.guide.startswith("Hamid"):
        #continue
    setName = row.guide.split("_")[0]
    realData[setName].append(int(row.modFreq))
    useScore = row.doench
    scoreData[setName].append((float(useScore), int(row.modFreq)))

    

print "Sets of guides:"
for setName, setVals in realData.items():
    print setName, setVals

#realData= {'snap25': [3, 2, 1], 'pTALRep': [1, 1, 3], 'pTAlReport1w': [3, 2, 1], 'cdk4int2': [3, 2, 1], 'cdk4int5': [3, 3, 1], 'Hamid': [2, 2, 3], 'SPARC': [3, 2, 2], 'OTR': [1, 1, 3]}

minSuccCount = 0
for setName, scoreList in scoreData.items():
    scoreList.sort()
    if scoreList[-1][1]==3:
        minSuccCount += 1
print "Total number of times the score is correct:", minSuccCount

totalSuccessCount = 0
n = 100000
for i in range(0, n):
    successCount = 0
    for gene, scores in realData.iteritems():
        random.shuffle(scores)
        score = scores[0]
        if score==3:
            successCount += 1
    if successCount >= minSuccCount:
        totalSuccessCount += 1

print totalSuccessCount, n
print "drawing a '3', >= 6 times, when drawing one from each set:", totalSuccessCount/float(n)
