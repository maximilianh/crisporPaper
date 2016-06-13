import glob
from collections import defaultdict, OrderedDict
import marshal

# - convert2.py creates the .../hart2016<cellLine><timepoint>.scores.tab files

def iterTsvRowsDict(ifh):
    " yield rows from a tab-sep table as OrderedDict "
    headers = ifh.readline().rstrip("\n").split("\t")
    for line in ifh:
        d = OrderedDict()
        row = line.rstrip("\n").split("\t")
        for name, val in zip(headers, row):
            d[name] = val
        yield d

dicts = list(iterTsvRowsDict(open("../hart2016Pseudo_hg19.scores.tab")))
seqToRow = {}
for d in dicts:
    seqToRow[d["seq"]] = d

modFreqs = marshal.load(open("modFreqs.marshal"))

for datasetName, seqChanges in modFreqs.iteritems():
    ofh = open("../%s.scores.tab" % datasetName, "w")
    headersDone = False
    geneCounts = defaultdict(int)
    for (gene, seq), foldChange in seqChanges.iteritems():
        geneCounts[gene]+=1
        d = seqToRow.get(seq)
        if d==None:
            continue
        d["modFreq"] = str(-float(foldChange))
        name = gene+"-"+str(geneCounts[gene])
        d["name"] = name
        d["dataset"] = datasetName
        row = d.values()
        row = [str(x) for x in row]
        if not headersDone:
            ofh.write("\t".join(d.keys())+"\n")
            headersDone = True
        ofh.write("\t".join(row))
        ofh.write("\n")
    print "wrote %s" % ofh.name

