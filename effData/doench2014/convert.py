from collections import defaultdict
import operator

from numpy import mean, log2, sum, array

import sys
sys.path.append("../../")
from annotateOffs import *

def parseExonIds(fname):
    " given a table seq,gene,exon,desc, return a dict with (gene, exonId) -> set of seqs (20mers)"
    ret = defaultdict(set)
    for row in iterTsvRows(fname):
        ret[(row.sym, row.exon)].add(row.seq[:20])
    return ret
        

geneExonSeqs = parseExonIds("../../out/seqToExon.tab")
# first get the raw abundances per 20mer guide sequence

lines = open("Doench-3026-S4.txt").read().splitlines()
heads = lines[1].split("\t")

# get the columns with the data
unsortCols = []
cd33Cols = []
cd13Cols = []
cd15Cols = []
for i, head in enumerate(heads):
    if head == "NB4 unsorted":
        unsortCols.append(i)
    if head == "NB4 CD33 neg":
        cd33Cols.append(i)
    if head == "MOLM13 CD15 neg":
        cd15Cols.append(i)
    if head == "NB4 CD13 neg":
        cd13Cols.append(i)

# convert the data to log abundances
guides = []
cd33Vals = []
cd13Vals = []
cd15Vals = []
allUnsortVals = []
seqs = []

for line in lines[2:]:
    fs = line.split("\t")
    seq = fs[0]
    trans = fs[1]
    #if trans!="ENST00000262262":
        #continue
    # get all unsorted and sorted values and their sums
    allUnsortVals.append( array([float(fs[col]) for col in unsortCols]) )
    cd33Vals.append( array([float(fs[col]) for col in cd33Cols]) )
    cd15Vals.append( array([float(fs[col]) for col in cd15Cols]) )
    cd13Vals.append( array([float(fs[col]) for col in cd13Cols]) )
    seqs.append((seq, trans))

cd33Sums = sum(cd33Vals, axis=0)
cd15Sums = sum(cd15Vals, axis=0)
cd13Sums = sum(cd13Vals, axis=0)
unsortSums = sum(allUnsortVals, axis=0)

rawCd33Scores = {}
rawCd15Scores = {}
rawCd13Scores = {}

guides = []

for (seq, trans), cd33Reads, cd15Reads, cd13Reads, unsortReads in zip(seqs, cd33Vals, cd15Vals, cd13Vals, allUnsortVals):
    cd33Logs = log2(1+1E6*(cd33Reads/cd33Sums))
    cd15Logs = log2(1+1E6*(cd15Reads/cd15Sums))
    cd13Logs = log2(1+1E6*(cd13Reads/cd13Sums))

    unsortLogs = log2(1+1E6*(unsortReads/unsortSums))

    cd33Score = mean(cd33Logs) - mean(unsortLogs)
    cd15Score = mean(cd15Logs) - mean(unsortLogs)
    cd13Score = mean(cd13Logs) - mean(unsortLogs)

    rawCd33Scores[seq] = 2**cd33Score
    rawCd15Scores[seq] = 2**cd15Score
    rawCd13Scores[seq] = 2**cd13Score
    # keep only CD33 guides for the ranking later
    if trans!="ENST00000262262":
        continue
    guides.append( (seq, 2**cd33Score) )
    #if seq=="TGGGGTGATTATGAGCACCG":
        #print "sums", sortSums, unsortSums
        #print unsortVals
        #print score

# add the percent-rank
#guides.sort(key=operator.itemgetter(1), reverse=True)
guides.sort(key=operator.itemgetter(1))
rankedGuides = []
guideData = {}
for i, (seq, val) in enumerate(guides):
    rankPerc = float(i)/len(guides)
    rankedGuides.append( (seq, val, rankPerc) )
    guideData[seq] = (val, rankPerc)

#ofh = open("doenchCD33.tab", "w")
#for seq, finalVal in guides:
    #ofh.write("%s\t%f\n" % (seq, finalVal))
#print 'wrote %s' % ofh.name

ofh = open("../doench2014-Hs.tab", "w")
ofh2 = open("../doench2014-Mm.tab", "w")
ofh3 = open("../doench2014-CD33Exon2.tab", "w")
ofh4 = open("../doench2014-CD33Exon3.tab", "w")
ofh5 = open("../doench2014-CD13Exon10.tab", "w")

headers = ["guide", "seq", "modFreq", "doenchSeq", "origDoenchScore", "CD33Abund", "CD15Abund", "CD13Abund", "geneAbund"]
ofh.write("\t".join(headers)+"\n")
ofh2.write("\t".join(headers)+"\n")

headers.append("origRankPerc")
headers.append("myRankPerc")
ofh3.write("\t".join(headers)+"\n")
ofh4.write("\t".join(headers)+"\n")
ofh5.write("\t".join(headers)+"\n")

#exon2Seqs = open("exon2Guides.txt").read().splitlines()
#print exon2Seqs
seqs3 = geneExonSeqs[("CD33", "exon2")]
seqs4 = geneExonSeqs[("CD33", "exon3")]
seqs5 = geneExonSeqs[("ANPEP", "exon10")] # ANPEP = CD13

geneCounts = defaultdict(int)
for line in open("Doench-nbt.3026-S5.txt"):
    if line.startswith("#"):
        continue
    fs = line.split()
    doenchSeq = fs[1]
    gene = fs[4]
    gene = gene.replace("H2-K", "H2K")
    # not "modification frequency" at all, but we just call it like this for now, 
    # so the file is compatible with the other files
    modFreq = fs[7]
    geneCounts[gene]+=1
    origDoenchScore = fs[-1]
    ensId = fs[5]
    if ensId.startswith("ENST"):
        rowOfh = ofh
    else:
        rowOfh = ofh2
    guideId = gene+"-"+str(geneCounts[gene])
    seq = doenchSeq[4:27]
    shortSeq = seq[:20]

    # output all three abundance scores, and also one column with the right one for this guide
    cd33 = rawCd33Scores.get(seq[:20], 0)
    cd15 = rawCd15Scores.get(seq[:20], 0)
    cd13 = rawCd13Scores.get(seq[:20], 0)

    rawAb = None
    if gene=="CD33":
        rawAb = cd33
    if gene=="CD13":
        rawAb = cd13
    if gene=="CD15":
        rawAb = cd15

    row = [guideId, seq, modFreq, doenchSeq, origDoenchScore, str(cd33), str(cd15), str(cd13), str(rawAb)]
    rowOfh.write("\t".join(row)+"\n")

    exOfh = None
    #print seqs3, seqs4, seqs5
    #print shortSeq
    if shortSeq in seqs3:
        exOfh = ofh3
    elif shortSeq in seqs4:
        exOfh = ofh4
    elif shortSeq in seqs5:
        exOfh = ofh5

    # the exon-specific files get raw abundance scores in the "modFreq" column
    if exOfh is not None:
        rawScore, myRankPerc = guideData.get(shortSeq, (None, None))
        row.append(row[2])
        row.append(str(myRankPerc))
        print "adding raw", row, rawAb
        row[2] = str(rawAb)
        exOfh.write("\t".join(row)+"\n")

ofh.close()
ofh2.close()

print "wrote output to %s and %s" % (ofh.name, ofh2.name)
print "wrote output to %s" % (ofh3.name)
print "wrote output to %s" % (ofh4.name)
print "wrote output to %s" % (ofh5.name)

#ofh = open("../doench2014-S10-414Genes.tab", "w")
#
#ofh.write("\t".join(headers)+"\n")
#
#geneCounts = defaultdict(int)
#for line in open("Doench-nbt.3026-S7.txt"):
#    if line.startswith("#"):
#        continue
#    fs = line.split()
#    doenchSeq = fs[1]
#    gene = fs[3]
#    gene = gene.replace("H2-K", "H2K")
#    # not "modification frequency" at all, but we just call it like this for now, 
#    # so the file is compatible with the other files
#    modFreq = fs[-2]
#    geneCounts[gene]+=1
#    origDoenchScore = fs[-1]
#    guideId = gene+"-"+str(geneCounts[gene])
#    seq = doenchSeq[4:27]
#    row = [guideId, seq, modFreq, doenchSeq, origDoenchScore]
#    ofh.write("\t".join(row)+"\n")
#ofh.close()
#
#print "wrote output to %s" % ofh.name

