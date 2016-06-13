# overwrite the mod frequencies from a .scores.tab file with the modification frequencies from another
# file

import sys
from collections import OrderedDict

def iterTsvRowsDict(ifh):
    " yield rows from a tab-sep table as OrderedDict "
    headers = ifh.readline().rstrip("\n").split("\t")
    for line in ifh:
        d = OrderedDict()
        row = line.rstrip("\n").split("\t")
        for name, val in zip(headers, row):
            d[name] = val
        yield d

scoresFname, modFreqFname, outFname = sys.argv[1:]

modFreqs = {}
for row in iterTsvRowsDict(open(modFreqFname)):
    modFreqs[row["guide"]] = row["modFreq"]
    
ofh = open(outFname, "w")
hDone = False
for row in iterTsvRowsDict(open(scoresFname)):
    if not hDone:
        ofh.write("\t".join(row.keys())+"\n")
        hDone = True

    newFreq = modFreqs.get(row["guide"])
    if newFreq==None:
        print "not found %s" % row["guide"]
        continue
    row["modFreq"] = newFreq
    ofh.write("\t".join(row.values())+"\n")

