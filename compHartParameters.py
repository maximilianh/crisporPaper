# obtain score <-> modFreq correlations for the hart2016 dataset using various filtering criteria
# a few different timepoints, various essential gene lists
# this was a try to figure out why the Hart correlations are sometimes better and sometimes worse against the
# the scores
import glob
import annotateOffs
from collections import defaultdict
from scipy.stats import pearsonr
import sys

seqToGene = {}

timeList = ["T12B", "T15B", "T18B"]

def parseTab(fname):
    " parse an R dataframe, return seq -> dict timepoint -> float "
    global seqToGene
    ifh = open(fname)
    print "reading", fname
    headers = ifh.readline().strip("\n").split("\t")
    tFields= [(x,y) for x, y in enumerate(headers) if y.startswith("T")]
    timepoints= [x for x in headers if x.startswith("T")]
    seqData = {}
    for line in ifh:
        fs = line.rstrip("\n").split("\t")
        gene, seq = fs[0].split("_")
        seqToGene[seq] = gene

        rowDict = {}
        for i, tName in tFields:
            rowDict[tName] = float(fs[i])
        seqData[seq] = rowDict
    return seqData, timepoints

def parseGenes(fname):
    # read the genes to filter on, return list
    syms = set()
    for line in open(fname):
        syms.add(line.strip().split()[0])
    return syms

def parseOnlyEffScores(inFname):
    """ return a dict seq -> scoreType -> score 
    """
    print "reading %s" % inFname

    scores = {}
    #freqs  = {}
    for row in annotateOffs.iterTsvRows(inFname):
        seq = row.seq
        scores[seq] = {}
        dataDict = row._asdict()
        for st in annotateOffs.scoreTypes:
            scores[seq][st] = float(dataDict[st])

        guideSeq = row.seq[:20]
        scores[seq]["finalGc6"] = int(annotateOffs.countFinalGc(guideSeq, 6)>=4)
        scores[seq]["finalGg"] = int(guideSeq[-2:]=="GG")
        scores[seq]["modFreq"] = float(row.modFreq)

        #freq = float(row.modFreq)
        #freqs[seq] = freq
        
    assert(len(scores)!=0)
    return scores

def addAvg(seqFoldChanges, timePoints):
    # add "avg" to seqFoldChanges
    for seq, timeToFold in seqFoldChanges.iteritems():
        s = 0
        for t in timePoints:
            s += timeToFold[t]
        timeToFold["avg"] = float(s)/len(timeList)

        s = 0
        for t in timePoints:
            if t.startswith("T18"):
                s += timeToFold[t]
        timeToFold["AvgLate"] = float(s)/len(timeList)
    return seqFoldChanges

def main():
    #geneSets = []
    #for fname in ["essentialSymsInFiveCellLines.txt", "core-essential-genes-sym_HGNCID", "essential_sym_hgnc.csv"]:
        #geneSet = parseGenes("effData/hart2016/essGenes/"+fname)
        #geneSets.append( (fname.split(".")[0]+"/"+str(len(geneSet)), geneSet) )

    rows = []
    for line in open("effData/hart2016/fileNames.txt"):
        if "DLD" in line:
            continue
        #if not "HeLa" in line:
            #continue
        foldChangeFname, datasetName, readFname = line.strip().split()
        print "dataset %s" % datasetName

        seqFoldChanges, timePoints1 = parseTab("effData/hart2016/TKOFoldChange/"+foldChangeFname)
        readCounts, timePoints2 = parseTab("effData/hart2016/readCounts/readcount-"+readFname)

        print "possible time points in folds: %s, time points in reads: %s" % (timePoints1, timePoints2)
        allTimePoints = set(timePoints1).intersection(timePoints2)

        # add "avg" to seqFoldChanges and readCounts
        seqFoldChanges = addAvg(seqFoldChanges, allTimePoints)
        readCounts = addAvg(readCounts, allTimePoints)

        print "all time points: ", allTimePoints
        allTimePoints.add("Avg")

        #headers = ["libraryName", "essentialGeneList", "minReads", "minFoldChange", "timepoint", "guideCountPass"]
        #headers = ["libraryName", "minReads", "minFoldChange", "timepoint", "guideCountPass"]
        headers = ["libraryName", "minReads", "timepoint", "guideCountPass"]

        #for geneSetName, genes in geneSets:
        for minReads in [0, 100, 300]:
            #for minFoldChange in [0.0, 5.0]:
            #for minFoldChange in [0.0, 5.0]:
                #for timepoint in ["T12B", "T15B", "T18B", "avg"]:
            for timepoint in allTimePoints:
                seqPredScores = parseOnlyEffScores("effData/%s%s.scores.tab" % (datasetName, timepoint))

                #row = [datasetName, geneSetName, minReads, minFoldChange, timepoint]
                row = [datasetName, minReads, timepoint]
                print "processing %s" % str(row)
                predScores = defaultdict(list)
                foldChanges = []
                #for seq, thisPredScores in seqFoldChanges.iteritems():
                for seq, thisPredScores in seqPredScores.iteritems():
                    # check gene
                    gene = seqToGene[seq]
                    #if gene not in genes:
                        #continue
                    # check foldChange
                    #foldChangeTimePoints = seqFoldChanges[seq]["modFreq"]
                    #if foldChangeTimePoints is None:
                        #print "%s has no fold change data but has been mapped" % seq
                        #continue

                    #foldChange = -foldChangeTimePoints[timepoint]
                    foldChange = seqPredScores[seq]["modFreq"]
                    #if foldChange < minFoldChange:
                        #continue
                    # check reads at T0
                    if int(readCounts[seq]["T0"]) < minReads:
                        continue

                    # seq is accepted:
                    # add all pred scores and the fold change
                    for scoreName, val in thisPredScores.iteritems():
                        predScores[scoreName].append(val)
                    foldChanges.append(foldChange)

                if len(predScores)==0:
                    #print "no data for %s" % str(row)
                    continue
                #assert(len(predScores)!=0)
                #assert(len(foldChanges)==len(predScores))
                # guideCountPass

                # some of the filtering critera leave us with 5 genes.
                # just skip these
                if len(foldChanges)<300:
                    continue
                row.append(len(foldChanges))

                for scoreName, predScoreList in predScores.iteritems():
                    #print "XX", scoreName, len(predScoreList)
                    assert(len(predScoreList)==len(foldChanges))
                    corr, pVal = pearsonr(predScoreList, foldChanges)
                    if scoreName not in headers:
                        headers.append(scoreName)
                    corr = "%0.3f" % corr
                    row.append(corr)
                row = [str(x) for x in row]
                rows.append(row)
                #print "\t".join(row)

    #ofh = open(sys.argv[1], "w")
    ofh = open("out/hartParams.tab", "w")
    ofh.write( "\t".join(headers))
    ofh.write( "\n")
    for row in rows:
        ofh.write( "\t".join(row))
        ofh.write( "\n")
    print "wrote to %s" % ofh.name

main()
