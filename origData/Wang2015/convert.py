import glob
from collections import defaultdict

# had wang data as pseudo-bed files copied from Fig1, convert them here to 
# a simpler one-table format

targets = []
for line in open("targets.txt"):
    targets.append(line.split()[0])

sums = defaultdict(int)
rows = []
for fname in glob.glob("*.bed"):
    for line in open(fname):
        fs = line.rstrip("\n").split()
        chrom, start, end, seq, score = fs[:5]
        fname = fname.split(".")[0]
        seqType = "off-target"
        if seq in targets:
            seqType = "on-target"
        score = float(score)
        row = [fname, seq.upper(), score, seqType]
        sums[fname] += score
        rows.append(row)


for row in rows:
    row[2] = str(row[2] / sums[row[0]])
    print "\t".join(row)
