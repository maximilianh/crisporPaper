# write all out of frame scores for data in effData/ to out/oofScores.tab

from annotateOffs import *
import os
import glob
from os.path import basename, splitext

resFh = open("out/oofScores.tab", "w")
resFh.write("seq\tlongSeq\toofScore\n")

for fname in glob.glob("effData/*.ext.tab"):
#for fname in glob.glob("effData/xu2015Train*.ext.tab"):
    print "processing ",fname
    if "schoenig" in fname: # cannot process a multispecies mixed dataset
        continue
    seqs = []
    ofh = open("/tmp/temp.bed", "w")
    modFreqs = []
    oofScores = []
    doneGuides = set()
    rowCount = 0
    inSeqs = set()
    for row in iterTsvRows(fname):
        if "position" not in row._fields:
            print "no genomic position for file %s" % fname
            assert(False)
        if row.guide in doneGuides:
            print "duplicated guide %s" % row.guide
            assert(False)
        if fname=="effData/doench2014-Mm.ext.tab":
            print "XX SKIPPING MOUSE DATA, too many duplicates"
            break
        chrom, startEnd, strand = row.position.split(":")
        start, end = startEnd.split("-")
        start = int(start)
        end = int(end)
        assert(end-start==34)
        if strand=="+":
            start = start - 16
            end = end + 30
        else:
            start = start - 30
            end = end + 16
        assert(end-start==80)

        if end==198711178 and "doench2014-Mm" in fname:
            # one sequence in Doench-Mm lies outside of chrom boundaries
            continue
        #if row.seq=="GGGAGGAGATAAGAAGAGAAAGG":
            #asdfdf
        inSeqs.add(row.seq)
        row = [ chrom, start, end, row.seq, "0", strand ]
        row = [str(x) for x in row]
        ofh.write("\t".join(row)+"\n")
        rowCount +=1

    ofh.close()
    print "Wrote %d rows to %s" % (rowCount, ofh.name)
    dataset = basename(fname).split(".")[0]
    db = datasetToGenome.get(dataset, "hg19")
    print "dataset %s, genome %s" % (dataset, db)
    cmd = "twoBitToFa -noMask -bed=/tmp/temp.bed /gbdb/%s/%s.2bit /tmp/temp.fa" % (db, db)
    assert(os.system(cmd)==0)

    count = 0
    for shortSeq, longSeq in parseFastaAsList(open("/tmp/temp.fa")):
        if shortSeq in inSeqs:
            inSeqs.remove(shortSeq)
        checkPos = longSeq.upper().index(shortSeq.upper())
        if checkPos!=20:
            print "one-bp shift in sequence %s, skipping it" % shortSeq
            continue
        #print checkPos, longSeq, shortSeq
        #assert(checkPos == 20)
        oofScore = calcMicroHomolScore(longSeq, 40)[1]
        #print longSeq[41:43]
        assert(longSeq[41]=="G" and longSeq[42]=="G")
        row = [shortSeq, longSeq, str(oofScore)]
        resFh.write("\t".join(row))
        resFh.write("\n")
        count += 1
    print "got %s long sequences" % count
    print inSeqs

print "wrote %s" % resFh.name
