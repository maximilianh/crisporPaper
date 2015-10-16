import sys
sys.path.append("..")
from annotateOffs import *

guideSeqs = {}
for line in open("suppTable1.txt"):
    #12 EMX1 GGCAGAAGCTGGAGGAGGAA GGG +
    fs = line.split()
    guideSeqs[(fs[1].upper(), fs[0])] = fs[2]

ofh = open("hsuSingle.tab", "w")

ofh.write('name\totName\tguideSeq\totSeq\tmismatches\treadFraction\tdiffLogo\n')
for line in open ("NIHMS539743-supplement-Supp_File_S2_singleMismatch.txt"):
    if line.startswith("#"):
        continue
    fs = line.split("\t")
    if "no transfection" in fs[0] or "Cas9 only" in fs[0]:
        continue
    guideId = fs[0].split("-")[0].upper()
    #print line, guideId
    guideSeq = guideSeqs[tuple(guideId.split("."))]
    freq = fs[7]
    if freq=="NaN":
        continue
    if "," in freq:
        freq = freq.split(",")[0].strip('"').strip("[")
    otSeq = fs[6]
    mmCount, diffLogo = countMms(guideSeq, otSeq)
    row = (guideId, fs[0], guideSeq, otSeq, str(mmCount), freq, diffLogo)
    ofh.write("\t".join(row)+"\n")

print "wrote %s" % ofh.name
