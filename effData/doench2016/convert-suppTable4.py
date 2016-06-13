from collections import defaultdict
import operator

from numpy import mean, log2, sum, array
import glob

for fname in glob.glob("suppTable4/*.txt"):
    if "Avana" in fname:
        continue
    print "infile %s" % fname
    base, lib, vec = fname.split(".")[0].split("-")
    outFname = "../doench2016-%s-%s_hg19.guides.tab" % (lib, vec)
    print "outfile %s" % outFname
    ofh = open(outFname, "w")
    ofh.write("guide\tseq\tmodFreq\n")
    guideId = 0
    skipCount = 0
    skipFcCount = 0
    for line in open(fname):
        if line.lower().startswith("vector") or line.lower().startswith("condition") or line.lower().startswith("construct") or line.lower().startswith("spacer") or line.lower().startswith("sequence"):
            continue
        if len(line)<3: # skip empty lines
            continue
        fields =  line.rstrip("\n").split("\t")
        if len(fields)==6:
            seq, pDna, rep1, rep2, rep3, rep4 = fields
            avg = (float(rep1)+float(rep2)+float(rep3)+float(rep4))/4
        elif len(fields)==4:
            seq, pDna, repA, repB = fields
            avg = (float(repA)+float(repB))/2
        else:
            print line
            assert(False)

        if avg==0.0:
            # pseudocount for 0 assay result
            avg = 0.001

        # ignore guides that were not there at all
        pDna = float(pDna)
        if pDna < 0.0000000001:
            skipCount += 1
            continue

        foldChange = avg / pDna

        #if foldChange == 0.0:
            #foldChange = 0.001
            #skipFcCount += 1
            #continue

        lfc = log2(foldChange)

        name = "%s%s%d" % (lib, vec, guideId)
        guideId += 1

        row = [name, seq, str(lfc)]
        ofh.write("\t".join(row)+"\n")

    print "not enough pDNA, skipped %d lines" % skipCount
    #print "exp. result is 0, skipped %d lines" % skipFcCount

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

