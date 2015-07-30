# take all of-targets and assign them to one exon

from annotateOffs import *
import os
import glob

ofh = open("/tmp/temp.bed", "w")

for fname in glob.glob("effData/*.ext.tab"):
    for row in iterTsvRows(fname):
        if "position" not in row._fields:
            break
        chrom, startEnd, strand = row.position.split(":")
        start, end = startEnd.split("-")
        row = [chrom, start, end, row.seq, "0", strand]
        ofh.write("\t".join(row)+"\n")
ofh.close()

ref2Sym = readDict("exonAnnot/ref2sym.tab")

selectFname = "exonAnnot/exons.bed"
inFname = ofh.name
outFname = "/tmp/temp.tab"
cmd = "overlapSelect -idOutput %s %s %s" % (selectFname, inFname, outFname)
assert(os.system(cmd)==0)

ofh2 = open("out/seqToExon.tab", "w")
ofh2.write("seq\tsym\texon\tothers\n")

for seq, exons in readDictList("/tmp/temp.tab").iteritems():
    symExons = set()
    sym = None
    exonId = None
    for ex in exons:
        refId, exonId = ex.split("|")
        sym = ref2Sym.get(refId, refId)
        symExons.add(sym+"|"+exonId)
    others = ",".join(symExons)
    row = [seq, sym, exonId, others]
    ofh2.write("\t".join(row)+"\n")

print "wrote %s" % ofh2.name
