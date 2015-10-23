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
    realData[row.guide.split("_")[0]].append(int(row.modFreq))

print realData

totalSuccessCount = 0
n = 10000
for i in range(0, n):
    successCount = 0
    for gene, scores in realData.iteritems():
        random.shuffle(scores)
        score = scores[0]
        if score==3:
            successCount += 1
    if successCount >= 6:
        totalSuccessCount += 1

print totalSuccessCount, n, totalSuccessCount/float(n)
