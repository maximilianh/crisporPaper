# change the modFreq field in ../varshney2015.scores.tab to the 2nd score field

newVals = {}
for line in open("Varshney2015_Supplemental_Table1.txt"):
    fs = line.strip().split()
    if len(fs)<4:
        continue
    name = fs[0]
    val = (fs[2], fs[3])
    newVals[name] = val

ofh = open("../varshney2015mutF0.scores.tab", "w")
for line in open("../varshney2015.scores.tab"):
    fs = line.strip().split()
    if fs[1] in newVals:
        fs[3]=newVals[fs[1]][0]
    else:
        print fs[1], "not found, keeping old val"

    ofh.write("\t".join(fs))
    ofh.write("\n")
print "wrote to %s" % ofh.name

ofh = open("../varshney2015mutF1.scores.tab", "w")
for line in open("../varshney2015.scores.tab"):
    fs = line.strip().split()
    if fs[1] in newVals:
        fs[3]=newVals[fs[1]][1]
    else:
        print fs[1], "not found, keeping old val"

    ofh.write("\t".join(fs))
    ofh.write("\n")
print "wrote to %s" % ofh.name

