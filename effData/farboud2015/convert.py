from collections import defaultdict

guideIds = defaultdict(int)
ofh = open("../farboud2015.tab", "w")
ofh.write("guide\tseq\tmodFreq\n")
for line in open("genetics.115.175166-1.txt"):
    if line.startswith("#EOF"):
        break
    if line.startswith("#"):
        continue
    fs = line.split()
    name = fs[0].replace("-", "")
    for i in range(len(fs)):
        val = fs[i]
        if len(val)==20:
            guideSeq = val
            pam = fs[i+1].strip("(").strip(")")
            freq = fs[i+3]
            float(freq) # make sure it's a number
    guideIds[name] += 1
    ofh.write("\t".join([name+"-"+str(guideIds[name]), guideSeq+pam, freq])+"\n")

print "wrote to %s" % ofh.name
