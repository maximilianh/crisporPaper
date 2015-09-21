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
columns = defaultdict(list)
for i, dataName in enumerate(heads[3:]):
    columns[dataName].append(i)

# convert the data to log abundances
seqs = [] # all seqs
readCounts = defaultdict(list) # dict dataset name -> lists of arrays, one per seq

# read all read counts as a list of arrays, one per dataset (=group of columns with same header)
for line in lines[2:]:
    fs = line.split("\t")
    seq = fs[0]
    trans = fs[1]
    #if trans!="ENST00000262262":
        #continue
    # get all unsorted and sorted values and their sums
    for datasetName, colList in columns.iteritems():
        readCounts[datasetName].append(array([float(fs[col+3]) for col in colList]) )
    seqs.append((seq, trans))

assert(len(seqs)==len(readCounts["NB4 CD33 neg"]))

# also get the sums of all columns, one per dataset
sums = {}
for datasetName, dataVals in readCounts.iteritems():
    sums[datasetName] = sum(dataVals, axis=0)

logMeans = defaultdict(dict) # sequence -> datasetName -> logScore

# now get the mean log2 score for each sequence, for all datasets
for datasetName, dataReads in readCounts.iteritems():
    dataSums = sums[datasetName]
    #print len(dataReads)
    #print len(dataSums)
    assert(len(seqs) == len(dataReads))
    for (seq, trans), seqReads in zip(seqs, dataReads):
        assert(len(seqReads)== len(dataSums))
        logScore = mean(log2(1+1E6*(seqReads/dataSums)))
        logMeans[seq][datasetName] = logScore

rawScores = defaultdict(dict) # cell type -> sequence -> logScore
for seq, dataMeans in logMeans.iteritems():
    for dataset in ['MOLM13 CD15 neg', 'NB4 CD33 neg', 'MOLM13 CD33 neg', 'TF1 CD13 neg', 'NB4 CD13 neg', 'TF1 CD33 neg']:
        dataStr = dataset.replace(" neg","").replace(" ","_")
        scoreDiff = dataMeans[dataset] - dataMeans[dataset.split()[0]+" unsorted"]
        rawScore = 2**scoreDiff
        rawScores[dataStr][seq] = rawScore

#guides.append( (seq, 2**cd33Score) )
##if seq=="TGGGGTGATTATGAGCACCG":
#    #print "sums", sortSums, unsortSums
#    #print unsortVals
#    #print score
#
## add the percent-rank
##guides.sort(key=operator.itemgetter(1), reverse=True)
#guides.sort(key=operator.itemgetter(1))
#rankedGuides = []
#guideData = {}
#for i, (seq, val) in enumerate(guides):
#    rankPerc = float(i)/len(guides)
#    rankedGuides.append( (seq, val, rankPerc) )
#    guideData[seq] = (val, rankPerc)
#
#ofh = open("doenchCD33.tab", "w")
#for seq, finalVal in guides:
    #ofh.write("%s\t%f\n" % (seq, finalVal))
#print 'wrote %s' % ofh.name

headers = ["guide", "seq", "modFreq", "doenchSeq", "origDoenchScore", "rankPercent"]

for dataset, seqScores in rawScores.iteritems():
    ofh = open("../doench2014-Hs-%s.tab" % dataset, "w")
    ofh.write("\t".join(headers)+"\n")
    #print seqScores

    geneCounts = defaultdict(int)
    for line in open("Doench-nbt.3026-S5.txt"):
        if line.startswith("#"):
            continue
        fs = line.split()
        doenchSeq = fs[1]
        gene = fs[4]
        gene = gene.replace("H2-K", "H2K")
        # so the file is compatible with the other files
        rankPercent = fs[7]
        geneCounts[gene]+=1
        origDoenchScore = fs[-1]
        ensId = fs[5]
        if not ensId.startswith("ENST"):
            continue
        guideId = gene+"-"+str(geneCounts[gene])
        seq = doenchSeq[4:27]
        shortSeq = seq[:20]
        rawScore = seqScores[shortSeq]

        if not gene in dataset:
            continue

        row = [guideId, seq, str(rawScore), doenchSeq, origDoenchScore, rankPercent]
        ofh.write("\t".join(row)+"\n")

    print "wrote output to %s" % (ofh.name)

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

