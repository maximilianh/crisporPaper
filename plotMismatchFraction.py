from annotateOffs import *

import numpy as np
import matplotlib.pyplot as plt

minFrac = 0.001
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

def plotBoxplots(fractions):
    " create boxplot of read fractions for each mismatch count "
    outfname = "mismatchFraction.pdf"
    #data = []
    xVals = defaultdict(list)
    yVals = defaultdict(list)
    labels = []
    maxMM = 7
    for mmCount in range(1,maxMM):
        otScores = fractions[mmCount]
        labels.append(str(mmCount)+" mismatches\n("+str(len(otScores))+" offtargets)")
        #data.append(otScores)
        for  name, otScore in otScores:
            study = name.split("_")[0]
            #xVals[study].append(mmCount-0.25+random.random()/2)
            xVals[study].append(mmCount)
            yVals[study].append(otScore)

    #plt.boxplot(data)
    colors = ["green", "blue", "black", "yellow", "red", "grey", "orange"]
    markers = ["o", "s", "+", ">", "<", "^", "x"]
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
        studyNames.append(study)
        i+=1
    #plt.xticks(range(1,maxMM), labels)
    plt.yticks(range(1,maxMM), labels)
    #plt.title("Off-target cleavage by number of mismatches")
    #plt.ylim((0,0.30))
    ax = plt.gca()
    plt.xlim((0.00,0.30))
    plt.ylim((0.5,6.5))
    for yGrid in range(1,6):
        ax.axhline(y=float(yGrid)+0.5, ls=":", linewidth=0.2, color="black")
    #plt.ylabel("Fraction of off-targets with indels")
    plt.xlabel("Indel frequency, when > 0.1%")

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


headers = ["name", "guideSeq", "otSeq", "guideGc", "readFraction", "mismatches", "otScore", "diffLogo"]

inFname = "offtargetsFilt.tsv"
#ignoreStudies = ["Hsu"]
ignoreStudies = []

targetSeqs = {}
inRows = []
otCounts = defaultdict(int)
for row in iterTsvRows(inFname):
    study = row.name.split("_")[0]
    if study in ignoreStudies:
        continue
    if row.type=="on-target":
        targetSeqs[row.name] = row.seq
    else:
        otCounts[row.name] += 1
        inRows.append(row)

rows = []
fractions = defaultdict(list)

datCount = 0
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

    #if not row.name.startswith("Tsai"):
    #   continue
    if float(row.score)>minFrac:
        fractions[mmCount].append((row.name, float(row.score)))
        datCount +=1

print "total points", datCount


rows.sort(key=operator.itemgetter(6))

ofh = open("annotFiltOfftargets.tsv", "w")
ofh.write( "\t".join(headers) )
ofh.write( "\n")
for row in rows:
    assert(len(row)==len(headers))
    row = [str(x) for x in row]
    ofh.write( "\t".join(row)+"\n")
print "wrote %s" % ofh.name
    
plotBoxplots(fractions)
