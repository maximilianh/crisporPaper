import os, logging
logging.basicConfig(loglevel=logging.INFO)
from os.path import isfile, splitext, join
from annotateOffs import *
from collections import defaultdict

logging.basicConfig(loglevel=logging.INFO)

from scipy.stats import linregress

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

def plotScores(ax, scores, scoreType, extFname, annotate, xLabel, doLegend=False, isXu=False, annotPos=None, diam=30):
    " create scatter plot "
    data = defaultdict(list)
    names = []
    regrX = []
    regrY = []
    print "REMOVING points where mod frequency is 0"
    for row in iterTsvRows(extFname):
        y = float(row.modFreq)
        if y==0.0:
            continue
        x = scores[row.extSeq][scoreType]
        names.append( (x, y, row.guide) )
        if not isXu or not row.guide.startswith("AAVS"):
            regrX.append(x)
            regrY.append(y)

        if isXu:
            if row.guide.startswith("AR"):
                data["AR"].append( (x, y) )
            elif row.guide.startswith("FOX"):
                data["FOXA1"].append( (x, y) )
            elif row.guide.startswith("AAVS"):
                data["AAVS"].append( (x, y) )
            else:
                assert(False)
        else:
            data["data"].append( (x,y) )

    if isXu:
        markers = {"AR" : "o", "FOXA1": "^", "AAVS":"x"}
    else:
        markers = {"data" : "o"}

    figs = []
    labels = []
    for title, tuples in data.iteritems():
        fig = ax.scatter(*zip(*tuples), \
            alpha=.5, \
            marker=markers[title], \
            s=diam)
        figs.append(fig)
        labels.append(title)

    if doLegend:
        ax.legend(figs,
               labels,
               scatterpoints=1,
               loc='upper left',
               ncol=2,
               fontsize=8)

    #fit = np.polyfit(regrX, regrY, 1)
    #m, b = fit
    #fit_fn = np.poly1d(fit) 
    ##plt.plot(x, y, 'yo', x, m*x+b, '--k') 
    #plt.plot(x,y, 'yo', x, fit_fn(x), '--k')
    #print "XX", m,b
    slope, intercept, r_value, p_value, std_err = linregress(regrX,regrY)
    print 'r^2 value', r_value**2
    print 'p_value', p_value
    print 'standard deviation', std_err
    line = slope*np.asarray(regrX)+intercept
    ax.plot(regrX,line, linestyle='--', color="orange")
    ax.annotate("r^2=%0.3f\npVa=%0.4f" % (r_value**2, p_value), xy=(0.75,0.05), fontsize=8, xycoords='axes fraction')
    #figs.append(studyFig)
    #plt.ylim(0,1.08)
    if annotate:
        for x, y, label in names:
           ax.annotate(
              label, fontsize=7, rotation=30, ha="right", rotation_mode="anchor",
              xy = (x, y), xytext = (0,0), alpha=0.9,
              textcoords = 'offset points', va = 'bottom')

    #plt.xlabel(scoreType+" score")
    ax.set_xlabel(xLabel)


def extendTabAddScores(extFname, scores, scoreNames, outFname):
    " add columns for efficiency scores to extFname "
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

def plotDataset(datasetName, ax, db, yLabel="Modification frequency", isXu=False, annotate=False, diam=30):
    inFname = join("effData/"+datasetName+".tab")
    extFname = extendTabAddContext(inFname, db)
    scores = addDoenchAndScs(extFname)

    extendTabAddScores(extFname, scores, ["ssc", "doench", "svm", "chariRaw"], "out/%s-compDoenchSsc.tsv" % datasetName)

    #plt.subplot(figRow,4,1)
    #ax[0].set_ylabel("Modification frequency")
    ax[0].set_ylabel(yLabel)
    ax[0].set_title(datasetName)
    plotScores(ax[0], scores, "svm", extFname, annotate, "SVM score from Wang et al. 2014", doLegend=isXu, isXu=isXu, annotPos=(0.2, 80), diam=diam)

    #plt.subplot(figRow, 4, 2)
    plotScores(ax[1], scores, "doench", extFname, annotate, "Score from Doench et al. 2014", isXu=isXu, annotPos=(-0.1, 95), diam=diam)

    #plt.subplot(figRow,4, 3)
    plotScores(ax[2], scores, "ssc", extFname, annotate, "SSC Score from Xu et al. 2015", isXu=isXu, annotPos=(-0.1, 90), diam=diam)

    #plt.subplot(figRow,4,4)
    plotScores(ax[3], scores, "chariRaw", extFname, annotate, "SVM Score from Chari et al. 2015", isXu=isXu, annotPos=(-0.1, 90), diam=diam)



def main():
    #"XuData/modFreq.tab" 
    #extendTabAddContext("temp.tab")
    plotFname = "out/compEffScores.pdf" 

    fig, axArr = plt.subplots(5, 4, sharey="row")
    fig.set_size_inches(18,20)

    #plotDataset("doench2014-Hs", axArr[4], "hg19")
    plotDataset("doench2014-Mm", axArr[4], "mm9")
    #plotDataset("xu2015", axArr[0], "hg19", isXu=True)
    #plotDataset("varshney2015", axArr[1], "danRer10")
    #plotDataset("gagnon2014", axArr[2], "danRer10")
    #plotDataset("xu2015Train", axArr[3], "danRer10", diam=1)
    #plotDataset("doench2014-S10-414Genes", axArr[4], "hg19", diam=1)
    #plotDataset("doench2014-S7-9Genes", 5, "hg19", annotate=False, diam=1)

    fig.tight_layout()
    fig.savefig(plotFname, format = 'pdf')
    fig.savefig(plotFname.replace(".pdf", ".png"))
    print "wrote plot to %s, added .png" % plotFname
    plt.close()

main()
