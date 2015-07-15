#plot for each guide: spec score with and without taking into account off-targets with >4 mismatches

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
import matplotlib.pyplot as plt
import numpy
import matplotlib.backends.backend_pdf as pltBack
from collections import defaultdict

from annotateOffs import *

ignoreStudies = ["Hsu", "Cho"] # ignore studies that only validated known off-targets

def parseSeqs(inFname):
    " parse guide seqs and their off-targets and index by name and mismatch count "
    # read guide sequences
    guideSeqs = {}
    for row in iterTsvRows(inFname):
        study = row.name.split("_")[0]
        if study in ignoreStudies:
            continue
        if row.type=="on-target":
            guideSeqs[row.name] = row.seq

    # now parse off-targets and index by name and mismatch count
    guideMms = defaultdict(dict)
    for row in iterTsvRows(inFname):
        if row.type=="off-target":
            if row.name not in guideSeqs:
                continue
            guideSeq = guideSeqs[row.name]
            mmCount, diffLogo = countMmsAndLogo(guideSeq, row.seq)
            guideMms[row.name].setdefault(mmCount, list()).append( (row.seq, diffLogo) )

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
            for otSeq, diffLogo in mmList:
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
    inFname = "offtargetsFilt.tsv"
    guideSeqs, guideMms = parseSeqs(inFname)

    rows = makeDataRows(guideSeqs, guideMms)

    outFname = 'out/specScoreComp.tsv'
    writeRows(rows, outFname)

    annotateXys = {} # dict with (x,y) -> name to annotate with labels
    figs = []
    for mmMax, column, color, marker in [(3, 1, "red", "o"), (4, 2, "black", "+"), (5, 3, "green", "^")]:
        xVals = []
        yVals = []
        for row in rows:
            # skip guides that don't have any data with 5 or 6 mismatches
            if row[0]=='Tsai_VEGFA_site2':
                annotateXys[(row[column], row[4])]=row[0]
            if row[-1]==row[-2]==0:
                continue
            xVals.append(row[column]) # first spec score
            yVals.append(row[4]) # 6MM spec score

        fig = plt.scatter(xVals, yVals, \
            alpha=.5, \
            color=color, \
            marker=marker, \
            s=30)
        figs.append(fig)

    plt.legend(figs,
           ["mismatches <= 3", "mismatches <= 4", "mismatches <= 5"],
           scatterpoints=1,
           loc='upper left',
           ncol=1,
           fontsize=10)

    for xy, guideName in annotateXys.iteritems():
       x, y = xy
       if float(x)>65 or y>25:
           plt.annotate(
              guideName, fontsize=9, rotation=30, ha="right", rotation_mode="anchor",
              xy = (x, y), xytext = (0,0), alpha=.5,
              textcoords = 'offset points', va = 'bottom')

    plt.xlabel("Spec. Score using off-targets with 3, 4 or 5 mismatches")
    plt.ylabel("Spec. Score using all off-targets (<= 6 mismatches)")
    outfname = "out/specScoreMMComp"
    plt.savefig(outfname+".pdf", format = 'pdf')
    plt.savefig(outfname+".png")
    print "wrote %s.pdf / .png" % outfname

main()
