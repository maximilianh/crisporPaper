from collections import defaultdict
import sys

# use validated freq data from supp table 5, adding sequences from other suppl files

#cellType = "Hap1"

#.ofh = open("convert.tab", "w")
ofh = sys.stdout

seqs = {}

def readSeqs():
    global seqs
    for guideName in ["HBB", "VEGFA"]:
        # mock file: get map name -> seq
        for line in open(guideName.lower()+"Mock.txt"):
            fs = line.rstrip("\n").strip().split()
            seq = fs[-1]
            name = fs[0]
            #if name.lower()=="on-target":
            name = guideName+"_"+name
            #print name
            assert(name not in seqs)
            seqs[name] = seq

def conv(guideName, cellType):
    # validated file: get map name -> frequency
    sumFreqs = 0.0
    rows = []
    for line in open(guideName.lower()+"Valid"+cellType+".txt"):
        fs = line.split()
        freq = float(fs[3].replace("%",""))/100
        name = fs[0]
        seqType = "off-target"
        if "On-target" in name:
            seqType = "on-target"
        sumFreqs+= freq
        row = ["Kim/"+cellType+"_"+guideName, seqs[guideName+"_"+name], freq, seqType]
        #row = ["Kim"+"_"+guideName, seqs[name], freq, seqType]
        rows.append(row)

    for row in rows:
        #row[2] = str(row[2] / sumFreqs)

        #row = [fname, seq, score, seqType]
        row = [str(x) for x in row]
        ofh.write( "\t".join(row))
        ofh.write("\n")

readSeqs()
conv("VEGFA", "Hap1")
conv("HBB", "Hap1")
conv("VEGFA", "K562")
conv("HBB", "K562")
#print ('wrote convert.tab')
