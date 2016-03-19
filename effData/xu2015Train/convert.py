# all scores are multiplied by -1 to be easier to compare with the other scores
# skipping all Mouse data

import glob
from math import *
from collections import defaultdict
ofh1 = open("../xu2015TrainHl60.ext.tab", "w")
ofh1.write("guide\tseq\textSeq\tmodFreq\tposition\n")
ofh2 = open("../xu2015TrainKbm7.ext.tab", "w")
ofh2.write("guide\tseq\textSeq\tmodFreq\tposition\n")

# the koike-yusa positions by Xu et al are incorrect
# We use the 20mer to lookup the original mm9 sequence from the Koike-Yusa table
# UAAARGH!!
#ifh = open("nbt.2800-S7.tsv")
# chr3	32712213	32712236	Actl6a	+	CAATGCCCTCCGCGTGCCCA	GGG	AATGCCCTCCGCGTGCCCA
#koikeYusaSeqs = {}
#for line in ifh:
    #chrom, start, end, gene, strand, seq, pam, revSeq = line.strip().split("\t")
    #if start=="Start":
        #continue
    #koikeYusaSeqs[ ( chrom, int(start)) ]= seq
    # koike yusa start is 1-based and indicates the start of the PAM
    #start = int(start)-1
    #if strand=="+":
        #start-=20
    #koikeYusaSeqs[ seq.upper() ]= (chrom, start, strand)

# now convert files
#ofh3 = open("../xu2015TrainMEsc.ext.tab", "w")
ofh3 = open("../xu2015TrainMEsc.guides.tab", "w")
ofh3.write("guide\tseq\tmodFreq\n")
#ofh3.write("guide\tseq\textSeq\tmodFreq\tposition\n")

geneCounts = defaultdict(int)
for fname in glob.glob("*.txt"):
    if fname=="log.txt" or fname=="README.txt":
        continue
    # ignore mouse guides
    #if "mEsc" in fname:
        #continue
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
        # use the average of the two cell line log values
        #logScore = -(s1+s2)/2
        if "mEsc" in fname:
            score1 = float(fs[5])
            score2 = float(fs[6])
            avgScore = log((exp(score1) + exp(score2))/2.0)
            ofh3.write("%s\t%s\t%s\n" % (guideId, shortSeq[:20], avgScore))
            #    chrom, start, strand = koikeYusaSeqs[ shortSeq[:20].upper() ]
            #    start = int(start)
            #    # need position of 34mer
            #    if strand=="+":
            #        start = origStart-4
            #        end   = origStart+30
            #    else:
            #        start = origStart-10
            #        end   = start+34

            #    #ofh3.write("%s\t%s\t%s\n" % (guideId, shortSeq[:20], avgScore))
            #    posStr = "%s:%d-%d:%s" % (chrom, start, end, strand)
            #   ofh3.write("%s\t%s\t%s\t%f\t%s\n" % (guideId, shortSeq, longSeq, hl60Score, posStr))
        else:
            hl60Score = float(fs[5])
            kbm7Score = float(fs[6])
            ofh1.write("%s\t%s\t%s\t%f\t%s\n" % (guideId, shortSeq, longSeq, hl60Score, posStr))
            ofh2.write("%s\t%s\t%s\t%f\t%s\n" % (guideId, shortSeq, longSeq, kbm7Score, posStr))
            #ofh2.write("%s\t%s\t%f\n" % (guideId, shortSeq, logScore))
ofh1.close()
ofh2.close()
print "output written to %s" % ofh1.name
print "output written to %s" % ofh2.name
print "output written to %s" % ofh3.name
#print "output written to %s" % ofh2.name
