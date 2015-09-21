seqs = {}
for line in open("chari-S6.txt"):
    if line.startswith("#"):
        continue
    fs = line.strip().split()
    seqId = fs[0].replace("_SPTargets","")
    seqs[seqId] = fs

for line in open("chari-S16.txt"):
    if line.startswith("#"):
        continue
    fs = line.strip().split()
    name = fs[0]
    otherFs = seqs[name]
    row = fs
    row.extend(otherFs)
    print "\t".join(row)


