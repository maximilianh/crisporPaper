# get all 23mers from the effData/ directory and write their extSeqs and positions to out/blatCache.tab
import glob
from annotateOffs import *
ofh = open("out/blatCache.tab", "w")
ofh.write("seq\textSeq\tposition\n")
for fname in glob.glob("effData/*.ext.tab"):
    for row in iterTsvRows(fname):
        if "position" in row._fields:
            pos = row.position
        else:
            continue
        seq, extSeq = row.seq, row.extSeq
        newRow = [seq, extSeq, pos]
        ofh.write("\t".join(newRow))
        ofh.write("\n")
print "wrote to %s" % ofh.name
