ofh = open("../schoenig.tab", "w")
ofh.write("guide\tseq\tmodFreq\n")

for line in open("schoenigData.txt"):
    if line.startswith("#"):
        continue
    fs = line.strip().split("\t")
    name = fs[0].replace(" ", "_")
    seq = fs[1]
    if len(seq)!=23:
        print "not 23bp long!", line
        continue
    modFreq = fs[3]
    if modFreq=="Na":
        continue
    ofh.write("\t".join([name,seq,modFreq]))
    ofh.write("\n")

print "wrote output to %s" % ofh.name

