import os, logging, operator
from os.path import isfile, splitext, join
from annotateOffs import *
from collections import defaultdict

from scipy.stats import linregress, binom

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
import matplotlib.pyplot as plt
import numpy as np

#logging.getLogger().setLevel(logging.DEBUG)

def compPairs(fname, minScoreDiff):
    print "Comparing guides from %s, minimum score difference %f" % (fname, minScoreDiff)
    byGene = defaultdict(list) # dict gene -> list of (guideName, modFreq, scores)
    for row in iterTsvRows(fname):
        if float(row.modFreq)==0.0:
            continue
        gene = row.guide.split("-")[0]
        scores = {}
        scores["doench"] = float(row.doench)
        scores["ssc"] = float(row.ssc)
        scores["svm"] = float(row.svm)
        chariRaw, chariRank = lookupchariScore(row.extSeq[4:27])
        scores["chariRaw"] = chariRaw
        scores["chariRank"] = chariRank
        byGene[gene].append( (row.guide, float(row.modFreq), scores) )

    # keep only genes with two guides
    twoGuides = dict()
    for gene, guideList in byGene.iteritems():
        if len(guideList)==2:
            twoGuides[gene]=guideList
        elif len(guideList)>2:
            guideList.sort(key=operator.itemgetter(1))
            twoGuides[gene]=(guideList[0], guideList[1])
        else:
            continue

    # for each gene, test if the order of the modFreq scores is the same as the order of the scores
    scoreNames = ["doench", "ssc", "svm", "chariRaw"]
    okCounts = defaultdict(int)
    for gene, guidePair in twoGuides.iteritems():
        guide1, guide2 = guidePair
        guide1Name, guide2Name = guide1[0], guide2[0]
        freq1, freq2 = guide1[1], guide2[1]
        if abs(freq2-freq1) < minScoreDiff:
            #print abs(freq2-freq1), "is <0.1"
            logging.debug("difference not high enough")
            continue
            
        okCounts["all"] += 1
        scores1, scores2 = guide1[2], guide2[2]
        logging.debug("guides (%s, %s), modFreq (%f, %f), doench (%f,%f), ssc (%f,%f)" % (guide1Name, guide2Name, freq1, freq2, scores1["doench"], scores2["doench"], scores1["ssc"], scores2["ssc"]))
        anyOk = False
        if freq2 > freq1:
            for scoreName in scoreNames:
                if scores2[scoreName] > scores1[scoreName]:
                    logging.debug( scoreName+ " OK")
                    okCounts[scoreName] += 1
                    anyOk = True
        else:
            for scoreName in scoreNames:
                if scores2[scoreName] < scores1[scoreName]:
                    logging.debug( scoreName+ " OK")
                    okCounts[scoreName] += 1
                    anyOk = True
        if not anyOk:
            logging.debug( "No score was OK")

    geneCount = okCounts["all"]
    print "total number of genes:", geneCount
    for scoreType, scoreCount in okCounts.iteritems():
        if scoreType=="all":
            continue
        pVal = binom.sf(scoreCount-1,geneCount,0.5)
        print "%s was correct %d times (p-Val %f)" % (scoreType, scoreCount, pVal)

def main():
    # parse the input file
    #compPairs("out/xu2015-compDoenchSsc.tsv", 20)
    compPairs("out/varshney2015-compEffData.tsv", 20)
    #compPairs("out/xu2015Train-compDoenchSsc.tsv", 0.9)
    #compPairs("out/doench2014-S7-9Genes-compDoenchSsc.tsv", 1.0)
    compPairs("out/doench2014-Hs-compEffData.tsv", 0.5)

main()
