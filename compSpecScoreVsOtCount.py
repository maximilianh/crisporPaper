# compare guide specificity scores against ot counts
# based on compareMitCrisporSpecScore.py
from annotateOffs import *
import matplotlib.pyplot as plt
from collections import defaultdict
from os.path import isfile
import pickle

# save time-intensive scores between invocations 
TMPFNAME = "/tmp/guideSpecScores.pickle"

# size expansion factor for bubbles
BUBBLEFAC = 200.0

def parseOtCounts(fname):
    " return a tuple of two dicts strongOtCount, weakOtCounts each guideName -> int "
    strongOffs = defaultdict(int)
    weakOffs = defaultdict(int)
    otShareSum = defaultdict(float)
    for row in iterTsvRows(fname):
        rf = float(row.readFraction)

        if rf>0.01:
            strongOffs[row.name]+=1
        #if rf>0.001:
        weakOffs[row.name]+=1

        otShareSum[row.name]+=rf
    return strongOffs, weakOffs, otShareSum

def main():
    maxMismatches = 4
    guideValidOts, guideSeqs = parseOfftargets("out/annotFiltOfftargets.tsv", maxMismatches, False, ["GG", "AG", "GA"])
    strongOtCounts, weakOtCounts, otShareSum = parseOtCounts("out/annotFiltOfftargets.tsv")

    if not isfile(TMPFNAME):
        crisporOffs = parseCrispor("crisporOfftargets", guideSeqs, 4)
        mitOffs = parseMit("mitOfftargets", guideSeqs)
        scoreCache = {}
    else:
        scoreCache = pickle.load(open(TMPFNAME))

    ofh = open("out/specScoreVsOtCount.tsv", "w")
    headers = ["guide", "CRISPORSpecScore", "MITSpecScore", "strongOtCount", "weakOtCount"]
    ofh.write("\t".join(headers)+"\n")
    xValsCrispor = []
    xValsMit = []
    yValsWeak = []
    yValsStrong = []
    areas = [] # size of the dots in the plot, one per xVal
    rows = []
    for guideName, guideSeq in guideSeqs.iteritems():
        if guideName in scoreCache:
            mitScore, crisporScore = scoreCache[guideName]
        else:
            mitScore = calcMitGuideScore_offs(guideSeq, mitOffs[guideSeq])
            crisporScore = calcMitGuideScore_offs(guideSeq, crisporOffs[guideSeq])
            scoreCache[guideName] = (mitScore, crisporScore)

        weakOtCount = weakOtCounts[guideName]
        strongOtCount = strongOtCounts[guideName]
        row = [guideName, crisporScore, mitScore, weakOtCount, strongOtCount]

        xValsCrispor.append(crisporScore)
        xValsMit.append(mitScore)
        yValsWeak.append(weakOtCount)
        yValsStrong.append(strongOtCount)
        areas.append(200.0*otShareSum[guideName])

        row = [str(x) for x in row]
        rows.append(row)

    rows.sort()

    for row in rows:
        ofh.write( "\t".join(row)+'\n')
    ofh.close()
    print "output written to %s" % ofh.name

    pickle.dump(scoreCache, open(TMPFNAME, "w"))

    figs = []
    weaks = plt.scatter(xValsCrispor, yValsWeak, \
        alpha=.5, \
        marker="o", \
        s=areas)

    #strongs = plt.scatter(xValsCrispor, yValsStrong, \
        #alpha=.5, \
        #marker="o", \
        #s=20)

    mitWeak = plt.scatter(xValsMit, yValsWeak, \
        alpha=.5, \
        marker="x", \
        s=areas)

    plt.xticks(range(0, 100, 10))

    legPlots = []
    for frac in [0.001, 0.005, 0.01, 0.05, 0.10, 0.30, 0.5, 0.7, 0.9]:
        legPlots.append(
            plt.scatter([],[], s=BUBBLEFAC*frac, edgecolors='none', marker="o"),

        )

    leg1 = plt.legend(legPlots, ["0.1%", "0.5%", "1%", "5%", "10%", "30%", "50%", "70%", "90%"],
           loc='upper right',
           ncol=1,
           fontsize=12, scatterpoints=1, title="Sum of\noff-target\nmodification\nfrequencies")

    plt.gca().add_artist(leg1)

    plt.legend([weaks, mitWeak],
           #["all off-targets", "off-targets <1%"],
           ["CRISPOR", "crispr.mit.org"],
           scatterpoints=1,
           loc='upper left',
           ncol=1,
           fontsize=12)
    plotFname = "out/specScoreVsOtCount.pdf"
    plt.xlabel("MIT Specificity Score of sgRNA")
    plt.ylabel("Number of off-targets of sgRNA")
    plt.ylim(0,40)
    plt.tight_layout()
    plt.savefig(plotFname, format = 'pdf')
    plt.savefig(plotFname.replace(".pdf", ".png"))
    print "wrote plot to %s, added .png" % plotFname
    plt.close()

main()
