import glob
from collections import defaultdict
import marshal

# export the hart data to various datasets: one for each cell line and timepoint
# - use a list of "essential genes"
# - add an additional timepoint, "avg", the difference between T0 and the average of all timepoints

# the way to create these files is as follows:
# - convert1.py creates a ../hart2016Pseudo.guides.tab file with just the sequences and names
# - then run ../effDataAddContext and ../effDataAddScores
# - convert2.py then creates the .../hart2016<cellLine><timepoint>.scores.tab files

def parseGenes(fname):
    # read the genes to filter on, return list
    syms = set()
    for line in open(fname):
        syms.add(line.strip().split()[0])
    return syms

def parseTab(fname, essGenes):
    " parse an R dataframe, return (gene, seq) -> dict timepoint -> float "
    #global seqToGene
    print "reading", fname
    ifh = open(fname)
    headers = ifh.readline().strip("\n").split("\t")
    tFields= [(x,y) for x, y in enumerate(headers) if y.startswith("T")]
    timepoints= [x for x in headers if x.startswith("T")]
    seqData = {}
    for line in ifh:
        fs = line.rstrip("\n").split("\t")
        gene, seq = fs[0].split("_")
        #seqToGene[seq] = gene

        rowDict = {}
        for i, tName in tFields:
            rowDict[tName] = float(fs[i])
        if gene not in essGenes:
            continue
        seqData[(gene, seq)] = rowDict
    return seqData, timepoints

def addAvg(seqFoldChanges, timePoints):
    # add "avg" to seqFoldChanges
    for seq, timeToFold in seqFoldChanges.iteritems():
        s = 0
        for t in timePoints:
            s += timeToFold[t]
        timeToFold["Avg"] = float(s)/len(timePoints)
    return seqFoldChanges

# the best parameters in compHartParams.py were:
# hart2016Hct1162lib2 essentialSymsInFiveCellLines/829     300      0.0           avg       1903      
essGenes = parseGenes("essGenes/essentialSymsInFiveCellLines.txt")

modFreqs = defaultdict(dict)
seqs = set()

for line in open("fileNames.txt"):
    inFname, datasetName, readCountName = line.strip().split()
    if "DLD" in line:
        continue

    foldFname = "TKOFoldChange/"+inFname
    seqFoldChanges, timePoints1 = parseTab(foldFname, essGenes)
    readCounts, timePoints2 = parseTab("readCounts/readcount-"+readCountName, essGenes)

    allTimePoints = set(timePoints1).intersection(timePoints2)

    # add "avg" to seqFoldChanges and readCounts
    seqFoldChanges = addAvg(seqFoldChanges, allTimePoints)
    readCounts = addAvg(readCounts, allTimePoints)

    print "possible time points in folds: %s, time points in reads: %s" % (timePoints1, timePoints2)
    allTimePoints = set(timePoints1).intersection(timePoints2)
    allTimePoints.add("Avg")

    seqIds = defaultdict(int)
    for timepoint in allTimePoints:
        #ofh = open("%s%s_hg19.guides.tab" % (datasetName, timepoint), "w")
        #ofh.write("guide\tseq\tmodFreq\n")

        for (gene, seq), foldChanges in seqFoldChanges.iteritems():
            #if int(readCounts[(gene, seq)][timepoint]) < 200:
                #continue
            seqs.add(seq)
            seqIds[gene]+=1
            name = gene+"-"+str(seqIds[gene])
            modFreqs[datasetName+timepoint][(gene, seq)]= foldChanges[timepoint]
            #row = [seq, name, str(foldChanges[timepoint])]
            #ofh.write("\t".join(row))

            #ofh.write("\n")
        #print "wrote %s" % ofh.name


        #i = 0
        #geneCounts = defaultdict(int)
        #for line in ifh:
            #fs = line.rstrip("\n").split("\t")
            #gene, seq = fs[0].split("_")
            ##val = fs[-1]
            #val = fs[firstTField]
            #geneCounts[gene] += 1
            #i += 1
            #ofh.write("%s-%d\t%s\t%s\n" % (gene, geneCounts[gene], seq, val))
            #print "wrote %s, %d guides, %d genes" % (ofh.name, i, len(geneCounts))
            #ofh.close()

ofh = open("modFreqs.marshal", "w")
modFreqs = dict(modFreqs)
marshal.dump(modFreqs, ofh)
print "wrote %s" % ofh.name

ofh = open("../hart2016Pseudo_hg19.guides.tab", "w")
ofh.write("guide\tseq\tmodFreq\n")
for s in seqs:
    ofh.write("%s\t%s\t1.0\n" % (s, s))

print "wrote %s" % ofh.name
