import sys
import glob
from collections import defaultdict

doneGuides = set()
for fname in glob.glob("specScores*.tab"):
    print "Parsing %s" % fname
    for line in open(fname):
        guide = line.split()[0][:20]
        assert(len(guide)==20)
        doneGuides.add(guide)

fname = sys.argv[1]

print "reading guides from allGuides.txt"
allGuides = set()
guideCounts = defaultdict(int)
for line in open(fname):
    guide = line.strip().split()[0]
    assert(len(guide)==20)
    allGuides.add(guide)
    guideCounts[guide[:20]]+=1

print "writing to jobFiles/"
gcCount =0
repCount = 0
nonUniqueCount = 0
chunkId = 0
doneCount = 0
i = 0
chunk = []
for guide in allGuides:
    assert(len(guide)==20)
    if guideCounts[guide[:20]] != 1:
        nonUniqueCount += 1
        continue
    if guide[:20].count("C") + guide[:20].count("G") >= 15:
        gcCount += 1
        continue
    # remove repeats
    if len(set(guide).intersection("actg"))!=0:
        repCount += 1
        continue
    if guide in doneGuides:
        doneCount += 1
        continue
    
    chunk.append(guide+"AGG")
    i += 1

    if len(chunk)==10:
        newFname = "jobFiles/%05d.fa" % chunkId
        ofh = open(newFname, "w")
        print "Writing %s" % newFname
        ofh.write(">chunk%d\n" % chunkId)
        ofh.write("NNNNNNN".join(chunk))
        ofh.write("\n")
        ofh.close()
        chunk = []
        chunkId += 1

print i, "sequences written"
print gcCount, "sequences removed due to high GC"
print repCount, "sequences removed due to repeats"
print nonUniqueCount, "sequences removed as they are non-unique"
print doneCount, "sequences removed as they were already done before"
