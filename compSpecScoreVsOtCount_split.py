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
            
def makePlot(ax, xVals, yVals, areas):

    (highX, highY, highZ), (lowX, lowY, lowZ) = splitXyzVals(xVals, yVals, areas, 0.005)

    #print lowX, lowY, lowZ, 200*lowZ
    #edgecolor='none', \
    ax.scatter(lowX, lowY, \
        alpha=.7, \
        marker="x", \
        edgecolor="black", \
        s=BUBBLEFAC*lowZ*SMALLFAC)

    #print highX, highY, highZ
    ax.scatter(highX, highY, \
        alpha=.7, \
        marker="o", \
        s=BUBBLEFAC*highZ)
    #plt.scatter(xVals, yVals, \
        #alpha=.7, \
        #marker="o", \
        #s=BUBBLEFAC*np.array(areas))

    ax.set_xticks(range(0, 101, 20))
    ax.set_xlim(0,100)
    ax.set_ylim(0,70)

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
            ax.scatter([],[], s=size, edgecolors='blue', marker=fracChar, lw=1),

        )

    # add legend
    leg1 = ax.legend(legPlots, ["0.1%", "0.5%", "1%", "5%", "10%", "30%", "50%", "70%", "90%"],
           loc='upper right',
           #bbox_to_anchor=(1.15, 1), loc=2, borderaxespad=0., \
           ncol=1,
           fontsize=10, scatterpoints=1, title="Sum of\noff-target\nmodification\nfrequencies")
    plt.setp(leg1.get_title(),fontsize='small')
    leg1.get_frame().set_linewidth(0.1)
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
    xVals = np.array(xVals)
    return xVals, yVals

def makeTwoSubplots(xValsMit, yValsWeak, areas, mitHistXVals, mitHistYVals, otCountsHistMit, suffix):

    f, axarr = plt.subplots(1, 2)

    # left subplot: scatter plot with bubble sizes
    f.set_size_inches(10,4)
    #f.subplots_adjust(top=0.10) # or whatever
    leg1 = makePlot(axarr[0], xValsMit, yValsWeak, areas)
    xlab = axarr[0].set_xlabel("%s Guide Specificity Score" % suffix)
    axarr[0].set_ylabel("Off-targets found per guide", color="black")
    axarr[0].set_xlim(0, 91)
    axarr[0].set_xticks(range(0, 91, 10))
    axarr[0].annotate('A', xy=(-.10, -.15), xycoords='axes fraction', fontsize=16,
                    horizontalalignment='right', verticalalignment='bottom')

    # right subplot: a histogram
    axarr[1].set_xlabel("%s Guide Specificity Score" % suffix)
    mitHistXVals = np.array(mitHistXVals)
    genomeBars = axarr[1].bar(mitHistXVals+2, mitHistYVals, 3, edgecolor='white', color="green" , lw=1)
    otBars = axarr[1].bar(mitHistXVals+3+2, otCountsHistMit, 3, edgecolor='white', color="blue" , lw=1)
    axarr[1].set_ylim(0,50)
    ylab = axarr[1].set_ylabel('Frequency', color="black")
    axarr[1].set_yticks(range(0, 51, 10))
    axarr[1].set_xticks(range(0, 101, 10))
    axarr[1].set_yticklabels(["%d%%" % x for x in range(0, 51, 10)])
    leg2 = axarr[1].legend( (otBars, genomeBars), ('Tested guides', 'All unique guides in human exons'), loc="upper left" )
    texts = leg2.get_texts()
    plt.setp(texts,fontsize='small')
    axarr[1].annotate('B', xy=(-.15, -.15), xycoords='axes fraction', fontsize=16,
                    horizontalalignment='right', verticalalignment='bottom')


    plt.tight_layout()
    plotFname = "out/specScoreVsOtCount-%s.pdf" % suffix
    #plt.savefig(plotFname, format = 'pdf')
    plt.savefig(plotFname, format = 'pdf', bbox_extra_artists=(leg1,xlab,ylab), bbox_inches='tight')
    plt.savefig(plotFname.replace(".pdf", ".png"), bbox_extra_artists=(leg1,), bbox_inches='tight')
    print "wrote plot to %s, added .png" % plotFname
    plt.close()
def main():
    maxMismatches = 4
    guideValidOts, guideSeqs = parseOfftargets("out/annotFiltOfftargets.tsv", maxMismatches, False, None)

    # get sum of off-target frequencies
    strongOtCounts, weakOtCounts, otShareSum = parseOtCounts("out/annotFiltOfftargets.tsv")

    crisporHistXVals, crisporHistYVals = parseSpecScores("wholeGenome/specScores.tab", "/tmp/crisporCache.pickle")
    mitHistXVals, mitHistYVals = parseSpecScores("seleniumMit/seqScores.txt", "/tmp/mitCache.pickle")
    assert(sum(crisporHistYVals)-100.0<0.01)
    assert(sum(mitHistYVals)-100.0<0.01)

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
    otCountsHistMit = [0] * 10
    otCountsHistCrispor = [0] * 10
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
        otCountsHistMit[mitScore/10] += 1
        otCountsHistCrispor[crisporScore/10] += 1
        areas.append(otShareSum[guideName])

        row = [str(x) for x in row]
        rows.append(row)

    rows.sort()

    # transform to frequencies in %
    otCountsHistMit = [100*x / float(len(guideSeqs)) for x in otCountsHistMit]
    assert(sum(otCountsHistMit)==100)

    otCountsHistCrispor = [100*x / float(len(guideSeqs)) for x in otCountsHistCrispor]
    assert(sum(otCountsHistCrispor)==100)

    for row in rows:
        ofh.write( "\t".join(row)+'\n')
    ofh.close()
    print "output written to %s" % ofh.name

    pickle.dump(scoreCache, open(TMPFNAME, "w"))

    makeTwoSubplots(xValsMit, yValsWeak, areas, mitHistXVals, mitHistYVals, otCountsHistMit, "MIT")

    makeTwoSubplots(xValsCrispor, yValsWeak, areas, crisporHistXVals, crisporHistYVals, otCountsHistCrispor, "CRISPOR")

main()
