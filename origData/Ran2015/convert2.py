import os
#ofh1 = open("sg1.temp.bed", "w")
#ofh2 = open("sg2.temp.bed", "w")
freqs = {}
headers = ["name", "seq", "score", "type"]
print "\t".join(headers)
for line in open("suppTable5.txt"):
    if not "SpCas9" in line:
        continue
    fs = line.split()
    pos = fs[3]
    freq = fs[-1]
    #print fs
    ns = fs[0].split("-")
    name = ns[0]+"-"+ns[1]
    #print "Ran_"+name, freq, pos
    chrom, startEnd = pos.split(":")
    start, end = startEnd.split("-")
    row = [chrom, int(start)-50, int(end)+50, fs[0]]
    row = [str(x) for x in row]
    #if "sg1" in line:
        #ofh = ofh1
    #elif "sg2" in line:
        #ofh = ofh2
    #else:
        #assert(False)
    freqs[fs[0]] = freq
    #ofh.write("\t".join(row))
    #ofh.write("\n")
    if "ontarget" in line:
        seqType = "on-target"
    else:
        seqType = "off-target"
    if freq=="N.D.":
        continue
    row = ["Ran_"+name, "noseq", freq, seqType]
    print "\t".join(row)

#ofh1.close()
#ofh2.close()

#cmd = "overlapSelect %s %s sg1.tmp -idOutput" % (ofh1.name, "Ran_EMX1-sg1.bed")
#assert(os.system(cmd)==0)

#cmd = "overlapSelect %s %s sg2.tmp -idOutput" % (ofh2.name, "Ran_EMX1-sg2.bed")
#assert(os.system(cmd)==0)

#for fname in ["sg1.tmp", "sg2.tmp"]:
    #for line in open(fname):
        #if line.startswith("#inId"):
            #continue
        #seq, name = line.strip().split()
        #print seq, name, freqs[name]
        
