from annotateOffs import *

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors

from os.path import isfile, join
import logging

#minFrac = 0.01

# an alternative version of this plot can be created by specifying the argument "supp"
doSupp = False

if len(sys.argv)>1 and sys.argv[1]=="supp":
    global doSupp
    doSupp = True

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

def plotFractions(fractions, mmAllCount, minFrac, baseOutName, otCount):
    """ create plot with read fractions for each mismatch count and total sites in genome from mmAllCount
    fractions:  mismatchCount -> list of (guideName, freq)
    """
    outfname = baseOutName+".pdf"
    # get the total OT count
    totalCount = 0
    maxMM = 7
    guideNames = set()
    for i in range(1, maxMM):
        print i, fractions[i]
        totalCount += len(fractions[i])
        for guideName, frac in fractions.iteritems():
            guideNames.add(guideName)
    print "check: totalCount of all validated off-targets= ", totalCount

    labels = []
    xVals = defaultdict(list)
    yVals = defaultdict(list)

    for mmCount in range(maxMM, 0, -1):
        otScores = fractions[mmCount]
        count = len(otScores)
        countPerc = (100.0*float(count)/totalCount)
        #avgHitCount = mmAllCount[mmCount] / len(guideNames)
        hitCount = mmAllCount[mmCount]
        label = "%s mismatches:    \n%d genome hits\n%d / %d offtargets (%0.1f %%)" % \
            (mmCount, hitCount, count, totalCount, countPerc)
        #labels.append(str(mmCount)+" mismatches: \n"+str(count)+" offtargets (%0.1f %%)" % countPerc)
        labels.append(label)

        for  name, otScore in otScores:
            study = name.split("_")[0]
            xVals[study].append(mmCount)
            yVals[study].append(otScore)

    colors = ["green", "orange", "black", "indigo", "red", "grey", "black", "black", "blue", "orange"]
    markers = ["o", "s", "v", ">", "<", "^", "+", "x", "o", "s"]
    studyNames = []

    studies = xVals.keys()
    studies.sort(reverse=True)

    i=0
    figs= []
    for study in studies:
        sXVals = xVals[study]
        xYVals = yVals[study]
        sXVals = [x - 0.3 + (0.1)*i for x in sXVals]
        sXVals, xYVals = xYVals, sXVals
        # linewidth=0 makes circles disappear
        edgecol = matplotlib.colors.ColorConverter().to_rgba(colors[i], alpha=0.5)
        fig = plt.scatter(sXVals, xYVals, alpha=0.5, marker=markers[i], color=colors[i], s=25, edgecolor=edgecol)
        figs.append(fig)
        #study = study.split("/")[0]
        studyNames.append(study)
        i+=1
    xLabels = ["%d%%" % int(100*x) for x in np.arange(0, 0.31, 0.05)]
    print xLabels
    plt.xticks(np.arange(0, 0.31, 0.05), xLabels)
    plt.yticks(range(maxMM,0,-1), labels)
    #plt.title("Off-target cleavage by number of mismatches")
    #plt.ylim((0,0.30))
    ax = plt.gca()
    #for tic in ax.xaxis.get_major_ticks():
        #tic.tick1On = tic.tick2On = False
    ax.yaxis.set_tick_params(width=0)
    if doSupp:
        plt.xlim((0.00,0.10))
    else:
        plt.xlim((0.00,0.30))
    plt.ylim((0.5,6.5))
    for yGrid in range(1,6):
        ax.axhline(y=float(yGrid)+0.58, ls=":", linewidth=0.2, color="black")
    #plt.ylabel("Fraction of off-targets with indels")
    label = "Modification frequency"
    #if minFrac!=0.0:
        #label += " > %0.2f%%" % (100*minFrac)
    plt.xlabel(label)

    plt.legend(reversed(figs),
           reversed(studyNames),
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

def indexOfftargets(inRows, minFrac, targetSeqs):
    " index offtargets by mismatchCount and return as dict mismatchCount -> (guideName, frequency) "
    fractions = defaultdict(list)
    datCount = 0
    for row in inRows:
        guideSeq = targetSeqs[row.name]
        otSeq = row.seq
        mmCount, diffLogo = countMms(guideSeq[:-3], otSeq[:-3])
        if float(row.score)>=minFrac:
            fractions[mmCount].append((row.name, float(row.score)))
            datCount +=1
    print "minimum frequency=%f: total off-targets %d" % (minFrac, datCount)
    return fractions, datCount

def parseRawOfftargets(inFname, onlyGuides = None):
    """ parse the raw list of off-targets, in the format of offtargets.tsv.
    returns list of rows and a dict guideName -> guideSeq 
    """
    targetSeqs = {}
    inRows = []
    for row in iterTsvRows(inFname):
        #if removeCellLine:
            # by removing the prefix before /, treat Kim's two cell lines as one experiment
            #study = row.name.split("_")[0].split("/")[0]
        if onlyGuides:
            if row.name not in onlyGuides:
                continue
        if row.type=="on-target":
            targetSeqs[row.name] = row.seq
        else:
            inRows.append(row)
    return inRows, targetSeqs

def main():
    inFname = "out/offtargetsFilt.tsv"
    onlyGuides = None
    if doSupp:
        inFname = "offtargets.tsv"
        onlyGuides = ["Tsai_HEK293_sgRNA4", "Tsai_VEGFA_site2"]

    inRows, targetSeqs = parseRawOfftargets(inFname, onlyGuides=onlyGuides)

    siteCountsByMismatch = readSiteCounts("crisporOfftargets", targetSeqs)
    fractions, otCount = indexOfftargets(inRows, 0.0, targetSeqs)
    plotFractions(fractions, siteCountsByMismatch, 0.0, "out/mismatchFraction-all", otCount)

    fractions, otCount = indexOfftargets(inRows, 0.01, targetSeqs)
    plotFractions(fractions, siteCountsByMismatch, 0.01, "out/mismatchFraction-min1Perc", otCount)

    fractions, otCount = indexOfftargets(inRows, 0.001, targetSeqs)
    plotFractions(fractions, siteCountsByMismatch, 0.001, "out/mismatchFraction-min01Perc", otCount)

main()
