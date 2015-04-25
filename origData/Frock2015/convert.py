def conv(fname, name, seqField, scoreField, filterSig=False):
    for i, line in enumerate(open(fname)):
        fs = line.split()
        if i==0:
            seqType = "on-target"
        else:
            seqType = "off-target"
        #print fs
        seq = fs[seqField]+fs[seqField+1][:3]
        #score = fs[6].replace("***")
        score = fs[scoreField]
        if filterSig and not "***" in score and seqType=="off-target":
            continue
        row = ["Frock_"+name, seq, score.replace("***",""), seqType]
        print "\t".join(row)

conv("suppTable3.txt", "RAG1A", 3, 6, filterSig=True)
conv("supplTable6.txt", "RAG1B", 4, 7)
conv("supplTable7a.txt", "EMX", 4, 7)
conv("supplTable7b.txt", "VEGFA", 4, 7)
