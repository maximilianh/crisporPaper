ofh1 = open("../concordet2-Hs.tab", "w")
ofh2 = open("../concordet2-Mm.tab", "w")
ofh3 = open("../concordet2-Rn.tab", "w")
ofh1.write("guide\tseq\tmodFreq\n")
ofh2.write("guide\tseq\tmodFreq\n")
ofh3.write("guide\tseq\tmodFreq\n")

# n10	hFSHD	5	24	1	6	AUA72	JBR3 p57	cTATAGGATCCACAGGGAGGgGG
for line in open("concordet2.txt"):
    if line.startswith("#"):
        continue
    fs = line.strip().split("\t")
    seq = fs[-1].upper()
    actLevel = fs[4]
    name = fs[1]+"_"+fs[2]
    row = [name, seq, str(int(actLevel)+1)]
    if(len(seq)!=23):
        print "not 23 bp long: %s" % line
        continue
    if fs[1].startswith("h"):
        ofh = ofh1
    elif fs[1].startswith("m"):
        ofh = ofh2
    elif fs[1].startswith("r"):
        ofh = ofh3
    ofh.write("\t".join(row)+"\n")

print "wrote %s" % ofh1.name
print "wrote %s" % ofh2.name
print "wrote %s" % ofh3.name
