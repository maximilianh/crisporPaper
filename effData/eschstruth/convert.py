ofh = open("../eschstruth.guides.tab", "w")
ofh.write("guide\tseq\tmodFreq\n")

for line in open("AlexisEschstruthCorrected.txt"):
    if line.startswith("#"):
        continue
    fs = line.strip().split("\t")
    name = fs[0]+"_"+fs[1]
    seq = fs[3][15:15+23]
    eff = fs[4]
    if eff == "++":
        modFreq = "3"
    elif eff == "+" or eff=="+-":
        modFreq = "2"
    elif eff == "-":
        modFreq = "1"
    else:
        assert(False)
    name = name.replace("'","P")
    row = [name, seq, modFreq]
    ofh.write("\t".join(row)+"\n")

print "wrote %s" % ofh.name
