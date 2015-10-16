# compare guide specificity scores against ot counts
# based on compareMitCrisporSpecScore.py
# XX currently recalculating the MIT scores. Maybe use the originals (once the site is working again)
from annotateOffs import *
import matplotlib.pyplot as plt
from collections import defaultdict
from os.path import isfile
import pickle
import numpy as np

# save time-intensive scores between invocations 
TMPFNAME = "/tmp/guideSpecScores.pickle"

# size expansion factor for bubbles
BUBBLEFAC = 200.0
# expansion factor for very small bubbles
SMALLFAC= 15.0

def parseOtCounts(fname):
    " return a tuple of three dicts strongOtCount, weakOtCounts, offtargetSum, each is : guideName -> float "
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

def splitXyzVals(xVals, yVals, zVals, zCutoff):
    """ split three lists into two sets of lists, depending on the yVal 
    stupid hack, numpy array would be one line.
    """
    x1, y1, z1 = [], [], []
    x2, y2, z2 = [], [], []
    for x,y,z in zip(xVals, yVals, zVals):
        if z <= zCutoff:
            x2.append(x)
            y2.append(y)
            z2.append(z)
        else:
            x1.append(x)
            y1.append(y)
            z1.append(z)
    return (x1, y1, np.array(z1)), (x2, y2, np.array(z2))
            
def makePlot(xVals, yVals, areas):

    (highX, highY, highZ), (lowX, lowY, lowZ) = splitXyzVals(xVals, yVals, areas, 0.005)

    #print lowX, lowY, lowZ, 200*lowZ
    #edgecolor='none', \
    plt.scatter(lowX, lowY, \
        alpha=.7, \
        marker="x", \
        edgecolor="black", \
        s=BUBBLEFAC*lowZ*SMALLFAC)

    #print highX, highY, highZ
    plt.scatter(highX, highY, \
        alpha=.7, \
        marker="o", \
        s=BUBBLEFAC*highZ)
    #plt.scatter(xVals, yVals, \
        #alpha=.7, \
        #marker="o", \
        #s=BUBBLEFAC*np.array(areas))

    plt.xticks(range(0, 101, 10))
    plt.xlim(0,100)
    plt.ylim(0,70)

    # invisible markers, needed for legend
    legPlots = []
    for frac in [0.005, 0.01, 0.05, 0.10, 0.30, 0.5, 0.7, 0.9]:
        fracChar = "o"
        size = BUBBLEFAC*frac
        if frac == 0.005:
            fracChar = "x"
            size *= SMALLFAC
        #if frac == 0.005:
            #fracChar = "x"
            #size *= 10

        legPlots.append(
            plt.scatter([],[], s=size, edgecolors='blue', marker=fracChar, lw=1),

        )

    # add legend
    leg1 = plt.legend(legPlots, ["0.1%", "0.5%", "1%", "5%", "10%", "30%", "50%", "70%", "90%"],
           #loc='upper right',
           bbox_to_anchor=(1.15, 1), loc=2, borderaxespad=0., \
           ncol=1,
           fontsize=10, scatterpoints=1, title="Sum of\noff-target\nmodification\nfrequencies")
    plt.setp(leg1.get_title(),fontsize='small')
    return leg1

def parseSpecScores(fname, cacheFname):
    """ parse a file with (seq,specScore) and return a list 0,10 with the percentage for each bin
    As this is somewhat slow, cache the result in a temp file.
    """
    print "parsing", fname
    if isfile(cacheFname):
        print "reading score histogram from temp file %s" % cacheFname
        return pickle.load(open(cacheFname))

    hist = [0] * 10
    totalCount = 0
    for line in open(fname):
        if "None" in line:
            continue
        score = int(line.rstrip("\n").split()[1])
        if score==100:
            score=99
        binIdx = score/10
        hist[binIdx]+=1
        totalCount += 1

    xVals = range(0, 100, 10)
    yVals = [100*(float(x)/totalCount) for x in hist]
    #print fname, xVals, yVals
    pickle.dump((xVals, yVals), open(cacheFname, "w"))
    return xVals, yVals

def main():
    maxMismatches = 4
    guideValidOts, guideSeqs = parseOfftargets("out/annotFiltOfftargets.tsv", maxMismatches, False, None)

    # get sum of off-target frequencies
    strongOtCounts, weakOtCounts, otShareSum = parseOtCounts("out/annotFiltOfftargets.tsv")

    histXVals, histYVals = parseSpecScores("wholeGenome/specScores.tab", "/tmp/crisporCache.pickle")
    mitHistXVals, mitHistYVals = parseSpecScores("seleniumMit/seqScores.txt", "/tmp/mitCache.pickle")

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

        weakOtCount   = weakOtCounts[guideName]
        strongOtCount = strongOtCounts[guideName]
        row = [guideName, crisporScore, mitScore, weakOtCount, strongOtCount]

        xValsCrispor.append(crisporScore)
        xValsMit.append(mitScore)
        yValsWeak.append(weakOtCount)
        yValsStrong.append(strongOtCount)
        areas.append(otShareSum[guideName])

        row = [str(x) for x in row]
        rows.append(row)

    rows.sort()

    for row in rows:
        ofh.write( "\t".join(row)+'\n')
    ofh.close()
    print "output written to %s" % ofh.name

    pickle.dump(scoreCache, open(TMPFNAME, "w"))

    plt.figure(figsize=(5,5))

    #axy1 = plt.subplot(121)
    leg1 = makePlot(xValsMit, yValsWeak, areas)
    xlab = plt.xlabel("MIT Specificity Score")
    plt.ylabel("Off-targets found per guide sequence", color="blue")

    ax1b = plt.twinx()
    ax1b.bar(mitHistXVals, mitHistYVals, 10, edgecolor='white', color="lightblue" , alpha=0.4, lw=1)
    ax1b.set_ylim(0,25)
    ylab = ax1b.set_ylabel('Frequency of specificity in exons (unique 20mers)', color="grey")

    #ax2 = plt.subplot(122, sharey=axy1)
    #plt.setp( ax2.get_yticklabels(), visible=False)
    #plt.ylim(0,60)

    plotFname = "out/specScoreVsOtCount-MIT.pdf"
    plt.savefig(plotFname, format = 'pdf', bbox_extra_artists=(leg1,xlab,ylab), bbox_inches='tight')
    plt.savefig(plotFname.replace(".pdf", ".png"), bbox_extra_artists=(leg1,), bbox_inches='tight')
    print "wrote plot to %s, added .png" % plotFname
    plt.close()

    # make the CRISPOR plot
    plt.figure(figsize=(5,5))
    leg1 = makePlot(xValsCrispor, yValsWeak, areas)

    xlab = plt.xlabel("CRISPOR Specificity Score")
    plt.ylabel("Off-targets found per guide sequence", color="blue")
    
    # add 2nd y axis and plot histogram
    ax2b = plt.twinx()
    ylab = ax2b.set_ylabel('Frequency of specificity in exons (unique 20mers)', color="grey")
    ax2b.set_ylim(0,20)
    ax2b.bar(histXVals, histYVals, 10, edgecolor='white', color="lightblue" , alpha=0.4, lw=1)

    #plt.tight_layout()

    #plt.tight_layout()
    # plt.subplots_adjust(hspace=0) # doesn't work

    plotFname = "out/specScoreVsOtCount-CRISPOR.pdf"
    plt.savefig(plotFname, format = 'pdf', bbox_extra_artists=(leg1,xlab,ylab), bbox_inches='tight')
    plt.savefig(plotFname.replace(".pdf", ".png"), bbox_extra_artists=(leg1,), bbox_inches='tight')
    print "wrote plot to %s, added .png" % plotFname
    plt.close()

main()
