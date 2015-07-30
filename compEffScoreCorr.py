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

def makePlot(xVals, yVals, areas, markerChar):
    weaks = plt.scatter(xVals, yVals, \
        alpha=.7, \
        edgecolor='none', \
        marker=markerChar, \
        s=areas)

    #strongs = plt.scatter(xValsCrispor, yValsStrong, \
        #alpha=.5, \
        #marker="o", \
        #s=20)

    #mitWeak = plt.scatter(xValsMit, yValsWeak, \
        #alpha=.5, \
        #marker="x", \
        #s=areas)

    plt.xticks(range(0, 101, 10))
    plt.xlim(0,100)
    plt.ylim(0,70)

    legPlots = []
    for frac in [0.001, 0.005, 0.01, 0.05, 0.10, 0.30, 0.5, 0.7, 0.9]:
        legPlots.append(
            plt.scatter([],[], s=BUBBLEFAC*frac, edgecolors='none', marker=markerChar),

        )

    #plt.gca().add_artist(leg1)

    #plt.legend([weaks, mitWeak],
           ##["all off-targets", "off-targets <1%"],
           #["CRISPOR", "crispr.mit.org"],
           #scatterpoints=1,
           #loc='upper left',
           #ncol=1,
           #fontsize=12)
    #plt.ylim(0,40)
    return legPlots

def parseSpecScores(fname):
    " parse a file with (seq,specScore) and return a list 0,10 with the counts for each bin "
    print "parsing", fname
    hist = [0] * 10
    totalCount = 0
    for line in open(fname):
        score = int(line.rstrip("\n").split()[1])
        if score==100:
            score=99
        binIdx = score/10
        hist[binIdx]+=1
        totalCount += 1

    xVals = range(0, 100, 10)
    yVals = [100*(float(x)/totalCount) for x in hist]
    return xVals, yVals

def main():
    maxMismatches = 4
    guideValidOts, guideSeqs = parseOfftargets("out/annotFiltOfftargets.tsv", maxMismatches, False, None)
    strongOtCounts, weakOtCounts, otShareSum = parseOtCounts("out/annotFiltOfftargets.tsv")

    histXVals, histYVals = parseSpecScores("wholeGenome/specScores.tab")

    if not isfile(TMPFNAME):
        crisporOffs = parseCrispor("crisporOfftargets", guideSeqs, maxMismatches)
        mitOffs = parseMit("mitOfftargets", guideSeqs)
        scoreCache = {}
    else:
        print "Not recalculating guide scores. Reading guide scores from %s" % TMPFNAME
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

    axy1 = plt.figure(figsize=(10,5))

    axy1 = plt.subplot(121)
    makePlot(xValsMit, yValsWeak, areas, "o")
    plt.xlabel("MIT Specificity Score")
    plt.ylabel("Number of off-targets")

    ax2 = plt.subplot(122, sharey=axy1)
    plt.setp( ax2.get_yticklabels(), visible=False)
    #plt.ylim(0,60)
    plt.ylabel("Off-targets found per guide sequence", color="blue")

    legPlots = makePlot(xValsCrispor, yValsWeak, areas, "o")

    # add legend
    leg1 = plt.legend(legPlots, ["0.1%", "0.5%", "1%", "5%", "10%", "30%", "50%", "70%", "90%"],
           #loc='upper right',
           bbox_to_anchor=(1.15, 1), loc=2, borderaxespad=0., \
           ncol=1,
           fontsize=10, scatterpoints=1, title="Sum of\noff-target\nmodification\nfrequencies")
    plt.setp(leg1.get_title(),fontsize='small')
    xlab = plt.xlabel("Specificity Score")
    
    # add 2nd y axis
    ax2 = plt.twinx()
    ax2.set_ylabel('Frequency of specificity in exons (unique 20mers)', color="grey")


    plt.bar(histXVals, histYVals, 10, edgecolor='white', color="lightblue" , alpha=0.5, lw=1)

    plt.tight_layout()

    #plt.tight_layout()
    # plt.subplots_adjust(hspace=0) # doesn't work

    plotFname = "out/specScoreVsOtCount.pdf"
    plt.savefig(plotFname, format = 'pdf', bbox_extra_artists=(leg1,xlab), bbox_inches='tight')
    plt.savefig(plotFname.replace(".pdf", ".png"), bbox_extra_artists=(leg1,), bbox_inches='tight')
    print "wrote plot to %s, added .png" % plotFname
    plt.close()

main()
