import glob
from collections import defaultdict

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

# the best parameters in compHartParams.py were:
# hart2016Hct1162lib2 essentialSymsInFiveCellLines/829     300      0.0           avg       1903      
for line in open("fileNames.txt"):
    inFname, datasetName, readCountName = line.strip().split()
    foldName = open("TKOFoldChange/"+inFname)
    #print "processing %s" % ifh.name

    foldChanges = parseTab(foldName)
    readCounts = parseTab("readCounts/readcount-"+readFname)
    # search for the first header that starts with T
    #headers = ifh.readline() # skip headers
    #headers = headers.rstrip("\n").split("\t")
    #firstTField = 0
    #for i, h in enumerate(headers):
        #if h.startswith("T"):
            #firstTField = i
            #break
    #assert(firstTField!=0)

    ofh = open("../%s_hg19.guides.tab" % datasetName, "w")
    ofh.write("guide\tseq\tmodFreq\n")
    i = 0
    geneCounts = defaultdict(int)
    for line in ifh:
        fs = line.rstrip("\n").split("\t")
        gene, seq = fs[0].split("_")
        #val = fs[-1]
        val = fs[firstTField]
        geneCounts[gene] += 1
        i += 1
        ofh.write("%s-%d\t%s\t%s\n" % (gene, geneCounts[gene], seq, val))
    print "wrote %s, %d guides, %d genes" % (ofh.name, i, len(geneCounts))
    ofh.close()
