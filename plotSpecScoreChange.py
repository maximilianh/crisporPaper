#plot for each guide: spec score change when changing mismatch count from 3 to 4 or 5

import matplotlib
import pickle
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
import matplotlib.pyplot as plt
import numpy
import matplotlib.backends.backend_pdf as pltBack
from collections import defaultdict

# save time-intensive score calculations between invocations 
TMPFNAME = "/tmp/specScoresByMms.pickle"

from annotateOffs import *

def parseSeqs(inFname):
    " parse guide seqs and their off-targets and index by name and mismatch count "
    # get guide sequences
    # get off-targets and index by name and mismatch count
    guideSeqs = {}
    guideMms = defaultdict(dict)
    for row in iterTsvRows(inFname):
        guideSeqs[row.name] = row.guideSeq
        mmCount = int(row.mismatches)
        guideMms[row.name].setdefault(mmCount, list()).append(row.otSeq)

    return guideSeqs, guideMms

def makeDataRows(guideSeqs, guideMms):
    """"
    return the rows for the tab-sep output file, format: name,specScore,lt4specScore,list of MMs
    """
    guideNames = sorted(guideSeqs.keys())

    rows = []
    for guideName in guideNames:
        row = [guideName, "", "", ""]
        # sum of hitscores: for <= 4MMs and all MMs
        lt3Sum, lt4Sum, lt5Sum, allSum = 0.0, 0.0, 0.0, 0.0

        mmInfo = guideMms[guideName]
        for mmCount in range(0, 7):
            mmList = mmInfo.get(mmCount, [])
            #row.append("%d:%d" % (mmCount, len(mmList)))
            row.append(len(mmList))

            guideSeq = guideSeqs[guideName]
            for otSeq in mmList:
                hitScore = calcHitScore(guideSeq[:20], otSeq[:20])

                allSum += hitScore
                if mmCount <=3:
                    lt3Sum += hitScore
                if mmCount <=4:
                    lt4Sum += hitScore
                if mmCount <=5:
                    lt5Sum += hitScore

        lt3SpecScore = calcMitGuideScore(lt3Sum)
        lt4SpecScore = calcMitGuideScore(lt4Sum)
        lt5SpecScore = calcMitGuideScore(lt5Sum)
        allSpecScore = calcMitGuideScore(allSum)

        row[1] = lt3SpecScore
        row[2] = lt4SpecScore
        row[3] = lt5SpecScore
        row[4] = allSpecScore
        rows.append(row)
    return rows

def writeRows(rows, outFname):
    ofh = open(outFname, 'w')
    headers = ["guideName", "specScoreUpToMM3", "specScoreUpToMM4", "specScoreUpToMM5", "specScoreUpToMM6", "MM0", "MM1", "MM2", "MM3", "MM4", "MM5", "MM6"]
    ofh.write("\t".join(headers)+"\n")
    for row in rows:
        row = [str(x) for x in row]
        ofh.write("\t".join(row))
        ofh.write("\n")
    ofh.close()
    print "wrote %s" % outFname

def main():
    inFname = "out/annotFiltOfftargets.tsv"
    guideSeqs, guideMms = parseSeqs(inFname)

    outFname = 'out/specScoreComp.tsv'
    rows = makeDataRows(guideSeqs, guideMms)
    writeRows(rows, outFname)

    scoreCache = defaultdict(dict)
    if isfile(TMPFNAME):
        scoreCache = pickle.load(open(TMPFNAME))
    else:
        crisporOffs = parseCrispor("crisporOfftargets", guideSeqs, 9999)

    annotateXys = []
    figs = []
    notShown = []

    doneSeqs = set()
    for guideName, guideSeq in guideSeqs.iteritems():
        if guideSeq in doneSeqs:
            continue
        doneSeqs.add(guideSeq)
        xVals, yVals = [], []
        #for maxMm in range(6,4, -1):
        for maxMm in range(6,3, -1):
            if maxMm in scoreCache[guideName]:
                specScore = scoreCache[guideName][maxMm]
            else:
                specScore = calcMitGuideScore_offs(guideSeq, crisporOffs[guideSeq], maxMm, 0.1, 1.0)
                scoreCache[guideName][maxMm] = specScore
            xVals.append(maxMm)
            yVals.append(specScore)
            if maxMm==4:
                annotateXys.append( (maxMm, specScore-0.5, guideName) )

        fig = plt.plot(xVals, yVals, \
            color="k", \
            marker="None")
        print guideName, xVals, yVals

    pickle.dump(scoreCache, open(TMPFNAME, "w"))

    print "Not shown, because no 5MM or 6MM values: %s" % (", ".join(set(notShown)))
    labels = ["3 mismatches", "4 mismatches", "5 mismatches", "6 mismatches"]
    xTicks = [3, 4, 5, 6]
    plt.xticks(xTicks, list(labels), rotation='vertical')
    plt.yticks(range(0,100,10))

    for x, y, guideName in annotateXys:
           plt.annotate(
              guideName, fontsize=9, ha="right", rotation_mode="anchor",
              xy = (x, y), xytext = (0,0),
              textcoords = 'offset points', va = 'bottom')

    #plt.xlabel("Spec. Score using off-targets with 3, 4 or 5 mismatches")
    plt.ylabel("Specificity Score")
    outfname = "out/specScoreMMComp"
    plt.xlim(3.5, 6.1)
    plt.ylim(0, 97)
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(5.5, 10)
    plt.tight_layout()
    plt.savefig(outfname+".pdf", format = 'pdf')
    plt.savefig(outfname+".png")
    print "wrote %s.pdf / .png" % outfname

main()
