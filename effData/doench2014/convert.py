from collections import defaultdict

ofh = open("../doench2014-S7-9Genes.tab", "w")

headers = ["guide", "seq", "modFreq", "doenchSeq", "origDoenchScore"]
ofh.write("\t".join(headers)+"\n")

geneCounts = defaultdict(int)
for line in open("Doench-nbt.3026-S5.txt"):
    if line.startswith("#"):
        continue
    fs = line.split()
    doenchSeq = fs[1]
    gene = fs[4]
    gene = gene.replace("H2-K", "H2K")
    # not "modification frequency" at all, but we just call it like this for now, 
    # so the file is compatible with the other files
    modFreq = fs[7]
    geneCounts[gene]+=1
    origDoenchScore = fs[-1]
    guideId = gene+"-"+str(geneCounts[gene])
    seq = doenchSeq[4:27]
    row = [guideId, seq, modFreq, doenchSeq, origDoenchScore]
    ofh.write("\t".join(row)+"\n")
ofh.close()

print "wrote output to %s" % ofh.name

ofh = open("../doench2014-S10-414Genes.tab", "w")

ofh.write("\t".join(headers)+"\n")

geneCounts = defaultdict(int)
for line in open("Doench-nbt.3026-S7.txt"):
    if line.startswith("#"):
        continue
    fs = line.split()
    doenchSeq = fs[1]
    gene = fs[3]
    gene = gene.replace("H2-K", "H2K")
    # not "modification frequency" at all, but we just call it like this for now, 
    # so the file is compatible with the other files
    modFreq = fs[-2]
    geneCounts[gene]+=1
    origDoenchScore = fs[-1]
    guideId = gene+"-"+str(geneCounts[gene])
    seq = doenchSeq[4:27]
    row = [guideId, seq, modFreq, doenchSeq, origDoenchScore]
    ofh.write("\t".join(row)+"\n")
ofh.close()

print "wrote output to %s" % ofh.name

