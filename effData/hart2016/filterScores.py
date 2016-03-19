# go over the files in mapped/, filter them by genes and write to ../
import glob
from os.path import basename

fileToReads = {}
for line in open("fileNames.txt"):
    inFname, datasetName, readCountName = line.strip().split()
    fileToReads[datasetName] = "readcount-"+readCountName

# read the genes to filter on
syms = set()
#for line in open("essGenes/essentialSymsInFiveCellLines.txt"):
#for line in open("essGenes/core-essential-genes-sym_HGNCID"):
for line in open("essGenes/essential_sym_hgnc.csv"):
    syms.add(line.strip().split()[0])
print "Got %d gene symbols to filter on" % len(syms)

for fname in glob.glob("mapped/*.scores.tab"):
    readCounts = {}
    ifh = open("readCounts/"+fileToReads[basename(fname).split(".")[0].split("_")[0]])
    headers = ifh.readline()
    for line in ifh:
        fs = line.rstrip("\n").split("\t")
        seq = fs[0].split("_")[1]
        val = int(fs[-1])
        readCounts[seq]=val

    ifh = open(fname)
    headerLine = ifh.readline()

    outFname = fname.replace("mapped/", "../")
    ofh = open(outFname, "w")
    ofh.write(headerLine)

    lcount = 0
    for line in ifh:
        fs = line.strip().split("\t")
        fs[3] = str(-float(fs[3]))
        geneSym = fs[1].split("-")[0]
        seq = fs[2]
        if geneSym in syms and readCounts[seq]>300:
            ofh.write("\t".join(fs))
            ofh.write("\n")
            lcount +=1
        
    print "wrote %d lines, file %s" % (lcount, ofh.name)
    ofh.close()

