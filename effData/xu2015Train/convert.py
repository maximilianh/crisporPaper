# all scores are multiplied by -1 to be easier to compare with the other scores
# skipping all Mouse data

import glob
from collections import defaultdict
ofh1 = open("../xu2015TrainHl60.ext.tab", "w")
#ofh2 = open("../xu2015Train.tab", "w")
ofh1.write("guide\tseq\textSeq\tmodFreq\tposition\n")
#ofh2.write("guide\tseq\tmodFreq\n")
ofh2 = open("../xu2015TrainKbm7.ext.tab", "w")
#ofh2 = open("../xu2015Train.tab", "w")
ofh2.write("guide\tseq\textSeq\tmodFreq\tposition\n")

geneCounts = defaultdict(int)
for fname in glob.glob("*.txt"):
    if fname=="log.txt":
        continue
    # ignore mouse guides
    if "mEsc" in fname:
        continue
    for line in open(fname):
        if line.startswith("#"):
            continue
        fs = line.rstrip("\n").split("\t")
        # ABT1    chr13   23423620        -       CCTTCACGTACGCAACCTGCTCAGCGCCTACGGCGAGGTG        -3.941506488    -3.835378137
        gene, chrom, origStart, strand = fs[:4]
        origStart = int(origStart)

        # need position of 34mer
        if strand=="+":
            start = origStart-4
            end   = origStart+30
        else:
            start = origStart-10
            end   = start+34

        assert(end-start==34)

        posStr = "%s:%d-%d:%s" % (chrom, start, end, strand)

        geneCounts[gene]+=1
        guideId = gene+"-"+str(geneCounts[gene])
        fullSeq = fs[4].upper()
        # strip first 6 basepairs to get a 34mer
        longSeq = fullSeq[6:]
        # from that, strip first 4 basepairs and last 7 to get a 23mer
        shortSeq = longSeq[4:-7]

        assert(len(longSeq)==34)
        assert(len(shortSeq)==23)
        assert(shortSeq.endswith("GG"))
        hl60Score = float(fs[5])
        kbm7Score = float(fs[6])
        # use the average of the two cell line log values
        #logScore = -(s1+s2)/2
        ofh1.write("%s\t%s\t%s\t%f\t%s\n" % (guideId, shortSeq, longSeq, hl60Score, posStr))
        ofh2.write("%s\t%s\t%s\t%f\t%s\n" % (guideId, shortSeq, longSeq, kbm7Score, posStr))
        #ofh2.write("%s\t%s\t%f\n" % (guideId, shortSeq, logScore))
ofh1.close()
ofh2.close()
print "output written to %s" % ofh1.name
print "output written to %s" % ofh2.name
#print "output written to %s" % ofh2.name
