# convert all effData/*.scores.tab files to effData/bed/*.bed files

import glob
from annotateOffs import *

for fname in glob.glob("effData/*.scores.tab"):

    outRows = []
    for row in iterTsvRows(fname):
        chrom, startEnd, strand = row.position.split(":")
        start, end  = startEnd.split("-")
        start = int(start)
        end = int(end)
        if strand=="+":
            guideStart = start+50-20
        else:
            guideStart = start+50
        guideEnd = guideStart+20

        row = [chrom, start, end, row.guide, "0", strand]
        outRows.append(row)

    outRows.sort()

    outFname = "effData/bed/"+basename(fname).split(".")[0]+".bed"
    ofh = open(outFname, "w")
    for row in outRows:
        row = [str(x) for x in row]
        ofh.write("\t".join(row))
        ofh.write("\n")
    ofh.close()
    print "wrote %s" % ofh.name
