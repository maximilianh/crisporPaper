from collections import defaultdict
import operator

from numpy import mean, log2

# first get the raw abundances per 20mer guide sequence

lines = open("Doench-3026-S4.txt").read().splitlines()
heads = lines[1].split("\t")

# get the columns with the data
unsortCols = []
sortCols = []
for i, head in enumerate(heads):
    if head == "NB4 unsorted":
        unsortCols.append(i)
    if head == "NB4 CD33 neg":
    #if head == "NB4 CD13 neg":
        sortCols.append(i)

# convert the data to log abundances
guides = []
allSortVals = []
allUnsortVals = []
seqs = []

for line in lines[2:]:
    fs = line.split("\t")
    seq = fs[0]
    trans = fs[1]
    if trans!="ENST00000262262":
        continue
    # get all unsorted and sorted values and their sums
    unsortVals.append( np.array([int(fs[col]) for col in unsortCols]) )
    sortVals.append( np.array([int(fs[col]) for col in sortCols]) )
    seqs.append(seq)
    #totalSum += sum(sortVals)

    #unsortMean = mean(unsortVals)
    #sortMean = mean(sortVals)
    #unsortSum = sum(unsortVals)
    #sortSum = sum(sortVals)
    ## transform them to log2 values
    #if sortSum==0.0:
    #    continue
    #totalSum = unsortSum+sortSum
    #unsortLogs = [log2(1+(1000000*float(val)/totalSum)) for val in unsortVals]
    ##unsortLogs = [log2(1+(1000000*float(val)/unsortSum)) for val in unsortVals]
    #sortLogs = [log2(1+(1000000*float(val)/totalSum)) for val in sortVals]
    ##sortLogs = [log2(1+(1000000*float(val)/sortSum)) for val in sortVals]
    #unsortMean = mean(unsortLogs)
    #sortMean = mean(sortLogs)
    #finalVal = sortMean - unsortMean
    #print "guide", seq
    #print "unsorted reads", unsortVals
    #print "sorted", sortVals
    #print "sums of unsorted/sorted/total", unsortSum, sortSum, totalSum
    #print "log2 of unsorted reads", unsortLogs
    #print "log2 of sorted reads", sortLogs
    #print "average of unsorted log2s", unsortMean
    #print "average of sorted log2s", sortMean
    #print "difference of averages", finalVal
    #print "2 ^ difference of averages", 2**finalVal
    #print "4 ^ difference of averages", 2**(2**finalVal)
    #print "sorted/unsorted", float(sortSum) / unsortSum
    #print "sorted avg/unsorted avg", mean(sortVals) / mean(unsortVals)
    #print "sorted sum/total sum", float(sum(sortVals)) / (sum(unsortVals)+sum(sortVals))
    ##print "logf", logf
    ##print "division of avg", 2**(sortMean/unsortMean)
    ##print "division of sum", sum(sortVals)/sum(unsortVals)
    ##print "final", mean([log2(x) for x in sortVals]) - mean([log2(x) for x in unsortVals])
    ##print "test", 2**(sortMean/unsortMean)
    #guides.append( [seq, sortVals] )

allSums = 
# add the percent-rank
#guides.sort(key=operator.itemgetter(1), reverse=True)
guides.sort(key=operator.itemgetter(1))
rankedGuides = []
guideData = {}
for i, (seq, val) in enumerate(guides):
    rankPerc = float(i)/len(guides)
    rankedGuides.append( (seq, val, rankPerc) )
    guideData[seq] = (val, rankPerc)

ofh = open("doenchCD33.tab", "w")
for seq, finalVal in guides:
    ofh.write("%s\t%f\n" % (seq, finalVal))
print 'wrote %s' % ofh.name

ofh = open("../doench2014-Hs.tab", "w")
ofh2 = open("../doench2014-Mm.tab", "w")
ofh3 = open("../doench2014-CD33Exon2.tab", "w")

headers = ["guide", "seq", "modFreq", "doenchSeq", "origDoenchScore", "myRawVal", "myRankPerc"]
ofh.write("\t".join(headers)+"\n")
ofh2.write("\t".join(headers)+"\n")

ofh3.write("\t".join(headers)+"\n")

exon2Seqs = open("exon2Guides.txt").read().splitlines()
#print exon2Seqs

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
    myVal, myRankPerc = guideData.get(shortSeq, ("0.0", "0.0"))
    row = [guideId, seq, modFreq, doenchSeq, origDoenchScore, str(myVal), str(myRankPerc)]
    rowOfh.write("\t".join(row)+"\n")

    if shortSeq in exon2Seqs:
        ofh3.write("\t".join(row)+"\n")

ofh.close()
ofh2.close()

print "wrote output to %s and %s" % (ofh.name, ofh2.name)
print "wrote output to %s" % (ofh3.name)

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

