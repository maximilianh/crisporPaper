# plot the quintiles of all scores against all dataset, like the figure in the Doench 2015 paper

import os, logging
logging.basicConfig(loglevel=logging.INFO)
from os.path import isfile, splitext, join
from annotateOffs import *
from collections import defaultdict

logging.basicConfig(loglevel=logging.INFO)

from scipy.stats import linregress
from scipy.stats import gaussian_kde
from numpy.random import normal
from numpy import arange

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
import matplotlib.pyplot as plt
import numpy as np

# sqlUcsc hgcentral -e 'select * from blatServers where db="mm9"'
blatServers = {
    "hg19": ("blat4a", "17779"),
    "danRer10" : ("blat4c", "17863"),
    "mm9" : ("blat4c", "17779")
}

def extend23Mers(seqs, db):
    """ extend 23mers to 30mers for eff score calculations, seqs is a dict id -> seq 
    return a dict seqId -> list of (seq, genomePositionString)
    """
    # write 23mers to fa file
    inNames = set()
    ofh = open("/tmp/temp.fa", "w")
    for seqId, seq in seqs.iteritems():
        ofh.write(">%s\n%s\n" % (seqId, seq))
        inNames.add(seqId)
    ofh.close()

    print "running BLAT, writing to /tmp/temp.bed"
    blatServer, port = blatServers[db]
    cmd = "gfClient %s.soe.ucsc.edu %s /gbdb/%s /tmp/temp.fa /tmp/temp.psl -minScore=20 -nohead -minIdentity=100 -maxIntron=0 -dots=1 ; pslToBed /tmp/temp.psl /tmp/temp.bed" % (blatServer, port, db)
    os.system(cmd)

    matches = defaultdict(list) # seqId -> list of (chrom, start, end, strand)
    for line in open("/tmp/temp.bed"):
        chrom, start, end, name, score, strand = line.split()[:6]
        if "_hap" in chrom or "random" in chrom:
            continue
        ##print int(end)-int(start)
        if (int(end)-int(start))!=23:
            #print "SKIP", line
            continue
        matches[name].append( (chrom, int(start), int(end), strand) )

    notFoundNames = inNames - set(matches)
    logging.warn("These sequences were not found with BLAT: %s" % ",".join(notFoundNames))
    #assert( len(seqs) == len(matches) )

    # write matches to temp.bed file
    # SSC needs extension by +7 bp of the end position
    # Doench needs extension -4 of the start and +3 of the end pos
    print "Creating /tmp/tempExt.bed with extended matches"
    ofh = open("/tmp/tempExt.bed", "w")
    positions = []
    for seqId, matchTuples in matches.iteritems():
        if len(matchTuples)>1:
            logging.error("Multiple matches for %s, will require manual selection" % seqId)
            logging.error("%s" % matchTuples)
        for matchTuple in matchTuples:
            chrom, start, end, strand = matchTuple
            if strand=="+":
                start = start - 4
                end = end + 7
            else:
                start = start - 7
                end = end + 4
            row = [chrom, str(start), str(end), seqId, "0", strand]
            ofh.write("\t".join(row)+"\n")
            positions.append( "%s:%d-%d:%s" % (chrom, start, end, strand))
    ofh.close()

    cmd = "twoBitToFa /gbdb/%s/%s.2bit -bed=/tmp/tempExt.bed /tmp/tempExt.fa" % (db, db)
    os.system(cmd)
    seqs = parseFastaAsList(open("/tmp/tempExt.fa"))
    assert(len(seqs)==len(positions))

    ret = defaultdict(list)
    for seqData, pos in zip(seqs, positions):
        seqId, seq = seqData
        ret[seqId].append( (seq, pos) )

    return ret

def extendTabAddContext(fname, db):
    """
    add a column to a tab file with a seq column: new column is called seqExt with 34 mers
    additional column is pos, with chrom:start-end:strand
    """
    newFname = fname.replace(".tab", ".ext.tab")
    if isfile(newFname):
        logging.info("not recreating %s, already exists. Delete to recreate" % newFname)
        return newFname

    seqs = dict()
    for row in iterTsvRows(fname):
        print "XX row", row
        seqs[row.guide] = row.seq
    seqPosDict = extend23Mers(seqs, db)

    ofh = open(newFname, "w")
    ofh.write("\t".join(row._fields)+"\t")
    ofh.write("extSeq\tposition\n")

    for row in iterTsvRows(fname):
        guideName = row.guide
        for seqPos in seqPosDict[guideName]:
            seq, pos = seqPos
            ofh.write("\t".join(row)+"\t"+seq+"\t"+pos+"\n")
    ofh.close()
    print "Wrote result to %s" % ofh.name
    return newFname

def addDoenchAndScs(fname):
    " given tab file with extSeq column, return scs and doench scores "
    seqs = []
    for row in iterTsvRows(fname):
        seqs.append(row.extSeq)
    return calcEffScores(seqs)

#binCount = 5
binCount = 4

def violin_plot(ax,data,pos, bp=False):
    '''
    create violin plots on an axis
    '''
    dist = max(pos)-min(pos)
    w = min(0.15*max(dist,1.0),0.5)
    for d,p in zip(data,pos):
        k = gaussian_kde(d) #calculates the kernel density
        m = k.dataset.min() #lower bound of violin
        M = k.dataset.max() #upper bound of violin
        x = arange(m,M,(M-m)/100.) # support for violin
        v = k.evaluate(x) #violin profile (density curve)
        v = v/v.max()*w #scaling the violin to the available space
        ax.fill_betweenx(x,p,v+p,facecolor='y',alpha=0.3)
        ax.fill_betweenx(x,p,-v+p,facecolor='y',alpha=0.3)
    if bp:
        ax.boxplot(data,notch=1,positions=pos,vert=1)

def emptyLists(c):
    " return a list of x empty lists "
    ret = []
    for i in range(0, c):
        ret.append([])
    return ret

        
def plotQuintiles(ax, scores, scoreType, extFname, title):
    " create barplots "
    xVals = []
    yVals = []
    for row in iterTsvRows(extFname):
        y = float(row.modFreq)
        x = scores[row.extSeq][scoreType]
        xVals.append(x)
        yVals.append(y)

    xVals = convToRankPerc(xVals)
    yVals = convToRankPerc(yVals)

    # bin based on x value
    scoreBins = emptyLists(binCount)
    binSize = 1.0/binCount
    for x, y in zip(xVals, yVals):
        scoreBin = int(x/binSize)
        print x, y, scoreBin
        scoreBins[scoreBin].append(y)

    #for i in range(0, binCount):
        #print i, len(scoreBins[i]), scoreBins[i][:3]

    binLabels = range(0, binCount)
    #violin_plot(ax, scoreBins, binLabels)

    # for each bin, bin again based on y-Value
    scoreActBins = [] # list of lists. One list per score bin. In this list, one list per activity bin
    for scoreBin, actVals in enumerate(scoreBins):
        actBins = emptyLists(binCount)
        for y in actVals:
            actBin = int(y / binSize) # should this be round() instead?
            actBins[actBin].append(y)

        # convert bins to frequencies
        actFreqs = []
        for actBin in actBins:
            if len(actVals)!=0:
                actFreq = len(actBin)/float(len(actVals))
            else:
                actFreq = 0.0
            actFreqs.append(actFreq)

        scoreActBins.append(actFreqs)
        print "actFreqs for bin", scoreBin, actFreqs

    colors = [0.9, 0.75, 0.5, 0.3, 0.1]
    colors = [str(x) for x in colors]
    currentSums = [0]*binCount
    for actBinIdx in reversed(range(0, binCount)):
        # get all vals for a given activity bin, plotting works by y-bin
        binVals = []
        for scoreBins in scoreActBins:
            binVals.append(scoreBins[actBinIdx])
        color = colors[actBinIdx]
        print "binVals for color", color, actBinIdx, binVals
        ax.bar(binLabels, binVals, color=color, linewidth=0, bottom=currentSums)
        currentSums = currentSums+np.array(binVals)

def extendTabAddScores(extFname, scores, scoreNames, outFname):
    " add columns for efficiency scores and write to extFname "
    #outFname = splitext(extFname)[0]+".scores.tab"
    ofh = None
    for row in iterTsvRows(extFname):
        if ofh==None:
            ofh = open(outFname, "w")
            ofh.write("\t".join(row._fields))
            ofh.write("\t")
            ofh.write("\t".join(scoreNames))
            ofh.write("\n")

        ofh.write("\t".join(row)+"\t")

        rowScores = []
        for name in scoreNames:
            rowScores.append(scores[row.extSeq][name])
        ofh.write("\t".join([str(x) for x in rowScores]))
        ofh.write("\n")
    ofh.close()
    print "wrote data to %s" % ofh.name

def plotDataset(datasetName, ax, db, title, addColLabels=False):
    inFname = join("effData/"+datasetName+".tab")
    extFname = extendTabAddContext(inFname, db)
    scores = addDoenchAndScs(extFname)

    extendTabAddScores(extFname, scores, ["ssc", "doench", "svm", "chariRaw"], "out/%s-compEffData.tsv" % datasetName)

    #ax[0].set_ylabel(yLabel)
    plotQuintiles(ax[0], scores, "svm", extFname, "SVM Score from Wang et al. 2014")

    plotQuintiles(ax[1], scores, "doench", extFname, "Score from Doench et al. 2014")
    plotQuintiles(ax[2], scores, "ssc", extFname, "SSC Score from Xu et al. 2015")

    plotQuintiles(ax[3], scores, "chariRaw", extFname, "SVM Score from Chari et al. 2015")

    if addColLabels:
        ax[0].set_title("SVM score from Wang et al. 2014")
        ax[1].set_title("Score from Doench et al. 2014")
        ax[2].set_title("Score from Xu et al. 2015")
        ax[3].set_title("Score from Chari et al. 2015")

    # put the row desc into the left border
    # http://stackoverflow.com/questions/25812255/row-and-column-headers-in-matplotlibs-subplots
    ax[0].annotate(title, xy=(0, 0.5), xytext=(-ax[0].yaxis.labelpad - 5, 0), \
       xycoords=ax[0].yaxis.label, textcoords='offset points', \
       #textcoords='offset points', \
       size='medium', ha='right', va='top')


def main():
    #"XuData/modFreq.tab" 
    #extendTabAddContext("temp.tab")
    plotFname = "out/effQuintiles.pdf" 

    fig, axArr = plt.subplots(7, 4, sharey="row", sharex="col")
    fig.set_size_inches(20,30)

    plotDataset("xu2015", axArr[0], "hg19", "Xu 2015, validation data\ntwo methods, FOX/AR sequencing, AAVS measured protein KO", addColLabels=True)
    plotDataset("varshney2015", axArr[1], "danRer10", "Varshney 2015\nZebrafish, special injection method?")
    plotDataset("gagnon2014", axArr[2], "danRer10", "Gagnon 2014, Zebrafish\nwhy so low?")
    plotDataset("xu2015Train", axArr[3], "danRer10", "Xu 2015 training data\ntaken from\nwhich other paper?")
    plotDataset("doench2014-S10-414Genes", axArr[4], "hg19", "Zhang 2014 (?)\n also used by Doench?")
    plotDataset("doench2014-Hs", axArr[5], "hg19", "Doench 2014, human")
    plotDataset("doench2014-Mm", axArr[6], "mm9", "Doench 2014, mouse")

    fig.tight_layout()
    fig.subplots_adjust(left=0.15, top=0.95)
    #fig.set_title(title)

    fig.savefig(plotFname, format = 'pdf')
    fig.savefig(plotFname.replace(".pdf", ".png"))
    print "wrote plot to %s, added .png" % plotFname
    plt.close()

if __name__=="__main__":
    main()
