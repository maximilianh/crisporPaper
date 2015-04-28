import glob

targets = []
for line in open("targets.txt"):
    targets.append(line.split()[1])

for fname in glob.glob("*.bed"):
    for line in open(fname):
        fs = line.rstrip("\n").split()
        chrom, start, end, seq, score = fs[:5]
        fname = fname.split(".")[0]
        seqType = "off-target"
        if seq in targets:
            seqType = "on-target"
        score = str(float(score)/100000)
        row = [fname.replace("cho", "Cho"), seq.upper(), score, seqType]
        print "\t".join(row)
