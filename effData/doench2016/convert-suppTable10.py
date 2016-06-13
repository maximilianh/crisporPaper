from collections import defaultdict
import operator

from numpy import mean, log2, sum, array
import glob

for fname in glob.glob("suppTable10/*Gecko2*.txt"):
    print "infile %s" % fname
    base, lib, vec = fname.split(".")[0].split("-")
    outFname = "../doench2016negSel-%s-%s_hg19.guides.tab" % (lib, vec)
    print "outfile %s" % outFname
    ofh = open(outFname, "w")
    ofh.write("guide\tseq\tmodFreq\n")
    guideId = 0

    for line in open(fname):
        if line.lower().startswith("vector") or line.lower().startswith("condition") or line.lower().startswith("construct") or line.lower().startswith("spacer") or line.lower().startswith("sequence"):
            #print line
            #sdfdf
            continue
        #if len(line)<3: # skip empty lines
            #continue
        fields =  line.rstrip("\n").split("\t")
        if len(fields)==6:
            seq, gene, rep1, rep2, rep3, rep4 = fields
            avgLfc = (float(rep1)+float(rep2)+float(rep3)+float(rep4))/4
        elif len(fields)==4:
            seq, gene, repA, repB = fields
            avgLfc = (float(repA)+float(repB))/2
        elif len(fields)==7:
            seq, gene, repHtA, repHtB, repHtC, repA, repB = fields
            avgLfc = (float(repHtA)+float(repHtB)+float(repHtC))/3
            #avgLfcA375 = (float(repA)+float(repB))/2
        else:
            print line
            print len(fields)
            print fields
            sdfdf

        name = "%s%s%d" % (lib, vec, guideId)
        guideId += 1

        row = [name, seq, str(avgLfc)]
        ofh.write("\t".join(row)+"\n")
