from annotateOffs import *

import numpy as np
import matplotlib.pyplot as plt

from os.path import isfile, join
import logging

#minFrac = 0.01

def countMms(string1, string2):
    " count mismatches between two strings "
    mmCount = 0
    string1 = string1.upper()
    string2 = string2.upper()
    diffLogo = []
    for pos in range(0, len(string1)):
        if string1[pos]!=string2[pos]:
            mmCount+=1
            diffLogo.append("*")
        else:
            diffLogo.append(".")
    return mmCount, "".join(diffLogo)

def readSiteCounts(crisporDir, guideNames):
    " return a dict with mismatchCount -> number of genome matches from dir "
    mmCounts = defaultdict(int)
    for guideName in guideNames:
        # special handling for Kim's two cell lines:
        # use only the Hap1 data, skip the K562 data
        # as they are identical
        guideName = guideName.replace("/Hap1", "")
        if "K562" in guideName:
            continue
        fname = join(crisporDir, guideName+".tsv")
        if not isfile(fname):
            logging.warn("Could not read %s -- YOU MUST ADD THIS FILE!" % fname)
            continue
        for line in open(fname):
            if line.startswith("guideId"):
                continue
            fs = line.split("\t")
            mm = int(fs[3])
            mmCounts[mm]+=1
    return mmCounts

def plotFractions(fractions, mmAllCount, minFrac, baseOutName):
    """ create plot with read fractions for each mismatch count and total sites in genome from mmAllCount
    """
    outfname = baseOutName+".pdf"
    # get the total OT count
    totalCount = 0
    maxMM = 7
    for i in range(1, maxMM):
        totalCount += len(fractions[i])
    print "check: totalCount of all validated off-targets= ", totalCount

    labels = []
    xVals = defaultdict(list)
    yVals = defaultdict(list)

    for mmCount in range(maxMM, 0, -1):
        otScores = fractions[mmCount]
        count = len(otScores)
        countPerc = (100.0*float(count)/totalCount)
        potCount = mmAllCount[mmCount]
        label = "%s mismatches\n%d genomic matches\n%d actual offtargets (%0.1f %%)" % \
            (mmCount, potCount, count, countPerc)
        #labels.append(str(mmCount)+" mismatches: \n"+str(count)+" offtargets (%0.1f %%)" % countPerc)
        labels.append(label)

        for  name, otScore in otScores:
            study = name.split("_")[0]
            xVals[study].append(mmCount)
            yVals[study].append(otScore)

    colors = ["green", "blue", "black", "yellow", "red", "grey", "orange", "violet"]
    markers = ["o", "s", "+", ">", "<", "^", "x", "+"]
    studyNames = []

    i=0
    figs= []
    for study, sXVals in xVals.iteritems():
        xYVals = yVals[study]
        sXVals = [x - 0.3 + (0.1)*i for x in sXVals]
        sXVals, xYVals = xYVals, sXVals
        # linewidth=0 makes circles disappear
        fig = plt.scatter(sXVals, xYVals, alpha=0.4, marker=markers[i], color=colors[i], s=20, edgecolor=colors[i])
        figs.append(fig)
        #study = study.split("/")[0]
        studyNames.append(study)
        i+=1
    #plt.xticks(range(1,maxMM), labels)
    plt.yticks(range(maxMM,0,-1), labels)
    #plt.title("Off-target cleavage by number of mismatches")
    #plt.ylim((0,0.30))
    ax = plt.gca()
    plt.xlim((0.00,0.30))
    plt.ylim((0.5,6.5))
    for yGrid in range(1,6):
        ax.axhline(y=float(yGrid)+0.5, ls=":", linewidth=0.2, color="black")
    #plt.ylabel("Fraction of off-targets with indels")
    label = "Modification frequency"
    if minFrac!=0.0:
        label += " > %0.2f%%" % (100*minFrac)
    plt.xlabel(label)

    plt.legend(figs,
           studyNames,
           scatterpoints=1,
           loc='upper right',
           ncol=3,
           fontsize=10)

    plt.tight_layout()
    plt.savefig(outfname)
    print "wrote %s" % outfname

    outfname = outfname.replace("pdf", "png")
    plt.savefig(outfname)
    print "wrote %s" % outfname
    plt.close()

def parseOfftargets(inFname, ignoreStudies):
    " parse the off-targets, count and return indexed "
    targetSeqs = {}
    inRows = []
    otCounts = defaultdict(int)
    for row in iterTsvRows(inFname):
        # by removing the prefix before /, treat Kim's two cell lines as one experiment
        study = row.name.split("_")[0].split("/")[0]
        if study in ignoreStudies:
            continue
        if row.type=="on-target":
            targetSeqs[row.name] = row.seq
        else:
            datasetName = row.name
            datasetName = datasetName.replace("/Hap1","").replace("/K562","")
            otCounts[datasetName] += 1
            inRows.append(row)
    return inRows, targetSeqs, otCounts

def makeOutRows(inRows, targetSeqs):
    rows = []
    guideNames = set()
    for row in inRows:
        guideSeq = targetSeqs[row.name]
        otSeq = row.seq
        mmCount, diffLogo = countMms(guideSeq[:-3], otSeq[:-3])
        if len(guideSeq)==23:
            otScore = calcHitScore(guideSeq[:-3], otSeq[:-3])
        else:
            otScore = "NA_not20mer"
        guideGc = gcCont(guideSeq)
        otRow = [row.name, guideSeq, otSeq, str(guideGc), row.score, str(mmCount), otScore, diffLogo]
        rows.append(otRow)
        guideNames.add(row.name)

    rows.sort(key=operator.itemgetter(6))
    return rows, guideNames

def indexOfftargets(inRows, minFrac, targetSeqs):
    " index offtargets by mismatchCount and return as dict mismatchCount -> (guideName, frequency) "
    fractions = defaultdict(list)
    datCount = 0
    for row in inRows:
        guideSeq = targetSeqs[row.name]
        otSeq = row.seq
        mmCount, diffLogo = countMms(guideSeq[:-3], otSeq[:-3])
        if float(row.score)>minFrac:
            fractions[mmCount].append((row.name, float(row.score)))
            datCount +=1
    print "minimum frequency=%f: total off-targets %d" % (minFrac, datCount)
    return fractions

def main():
    headers = ["name", "guideSeq", "otSeq", "guideGc", "readFraction", "mismatches", "otScore", "diffLogo"]
    inFname = "offtargetsFilt.tsv"
    ignoreStudies = []

    inRows, targetSeqs, otCounts = parseOfftargets(inFname, ignoreStudies)

    rows, guideNames = makeOutRows(inRows, targetSeqs)

    # write out rows
    ofh = open("annotFiltOfftargets.tsv", "w")
    ofh.write( "\t".join(headers) )
    ofh.write( "\n")
    for row in rows:
        assert(len(row)==len(headers))
        row = [str(x) for x in row]
        ofh.write( "\t".join(row)+"\n")
    print "wrote %s" % ofh.name

    siteCountsByMismatch = readSiteCounts("crisporOfftargets", guideNames)
    fractions = indexOfftargets(inRows, 0.0, targetSeqs)
    plotFractions(fractions, siteCountsByMismatch, 0.0, "mismatchFraction-all")

    fractions = indexOfftargets(inRows, 0.01, targetSeqs)
    plotFractions(fractions, siteCountsByMismatch, 0.01, "mismatchFraction-min1Perc")

main()
