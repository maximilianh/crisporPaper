import glob
from collections import defaultdict
ofh = open("../xu2015Train.ext.tab", "w")
#ofh2 = open("../xu2015Train.tab", "w")
ofh.write("guide\tseq\textSeq\tmodFreq\n")
#ofh2.write("guide\tseq\tmodFreq\n")

geneCounts = defaultdict(int)
for fname in glob.glob("*.txt"):
    if fname=="log.txt":
        continue
    for line in open(fname):
        if line.startswith("#"):
            continue
        fs = line.rstrip("\n").split("\t")
        # ABT1    chr13   23423620        -       CCTTCACGTACGCAACCTGCTCAGCGCCTACGGCGAGGTG        -3.941506488    -3.835378137
        gene = fs[0]
        geneCounts[gene]+=1
        guideId = gene+"-"+str(geneCounts[gene])
        fullSeq = fs[4].upper()
        # strip first 7 basepairs
        longSeq = fullSeq[6:]
        # strip first 4 basepairs and last 7
        shortSeq = longSeq[4:-7]

        assert(len(longSeq)==34)
        print shortSeq, len(shortSeq)
        assert(len(shortSeq)==23)
        assert(shortSeq.endswith("GG"))
        s1 = float(fs[5])
        s2 = float(fs[6])
        logScore = (s1+s2)/2
        ofh.write("%s\t%s\t%s\t%f\n" % (guideId, shortSeq, longSeq, logScore))
        #ofh2.write("%s\t%s\t%f\n" % (guideId, shortSeq, logScore))
ofh.close()
print "output written to %s" % ofh.name
#print "output written to %s" % ofh2.name


        
