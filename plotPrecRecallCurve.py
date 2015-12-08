# for each dataset, get the top25% and the bottom25% and assign labels 1.0 and 0.0
# then get all scores for these sequences and write scores and labels to a file precRec/dataset-scoreType.tab
# then use the R ROCR package and write the x and y values to a file
# precRec/dataset-scoreType.precRec.tab 
# finally plot all these curves with matplotlib

import os, logging
from annotateOffs import *
from collections import defaultdict

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
    scoreNames = seqScores.values()[0]

    nameToScores = {}
    for scoreName in scoreNames:
        scores = []
        for s in seqs:
            scores.append(seqScores[s][scoreName])
        nameToScores[scoreName] = scores

    return nameToScores

seqScores, freqs = parseEffScores("chari2015Train")
topSeqs, bottomSeqs = getQuartiles(freqs)

labels = [1]*len(topSeqs)
labels.extend( [0]*len(bottomSeqs) )
allSeqs = topSeqs
allSeqs.extend(bottomSeqs)
allScores = getScores(seqScores, allSeqs)

print topScores



