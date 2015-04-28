from collections import defaultdict

def conv(fname, name, seqField, scoreField, filterSig=False):
    rows = []
    sums = defaultdict(float)
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
        score = score.replace("***", "")
        row = [name, seq, score, seqType]
        rows.append(row)
        if score!="NA":
            sums[name]+=float(score)

    sums = dict(sums)
    for row in rows:
        #print row
        #print sums
        if row[2]!="NA":
            row[2] = str(float(row[2])/sums[row[0]])
        row[0] = "Frock_"+row[0]
        print "\t".join(row)

# remove non-sign interactions?
conv("suppTable3.txt", "RAG1A", 3, 6, filterSig=False)
conv("supplTable6.txt", "RAG1B", 4, 7)
conv("supplTable7a.txt", "EMX", 4, 7)
conv("supplTable7b.txt", "VEGFA", 4, 7)
