seqs = {}
for line in open("chari-S6.txt"):
    if line.startswith("#"):
        continue
    fs = line.strip().split()
    seqId = fs[0].replace("_SPTargets","")
    seqs[seqId] = fs[1]

ofh = open("../chari2015Train293T.tab", "w")
ofh.write("guide\tseq\tmodFreq\n")
#for line in open("chari-S16.txt"):
for line in open("chari2015-suppData10-S13.txt"):
    if line.startswith("#"):
        continue
    fs = line.strip().split()
    modFreq = float(fs[2]) - float(fs[1])
    #modFreq = float(fs[2])
    if modFreq < 0.0:
        modFreq = 0.0
    name = fs[0]
    ofh.write("%s\t%s\t%s\n" % (name, seqs[name], str(modFreq)))

print "wrote to %s" % ofh.name

ofh = open("../chari2015TrainK562.tab", "w")
ofh.write("guide\tseq\tmodFreq\n")
for line in open("chari-S16.txt"):
    if line.startswith("#"):
        continue
    fs = line.strip().split()
    modFreq = float(fs[2]) - float(fs[1])
    #modFreq = float(fs[2])
    if modFreq < 0.0:
        modFreq = 0.0
    name = fs[0]
    ofh.write("%s\t%s\t%s\n" % (name, seqs[name], str(modFreq)))

print "wrote to %s" % ofh.name
