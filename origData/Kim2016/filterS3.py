# filter table S3 and keep only the Tsai sequences from it.

tsaiSeqs = open("tsaiSeqs.txt").read().splitlines()

for line in open("kim2016-suppTable3-hela.tsv"):
    if line.startswith("#"):
        print line,
    seq = line.split()[3]
    if seq.upper() not in tsaiSeqs:
        continue

    mmCount = 0
    for c in seq:
        if c.islower():
            mmCount +=1

    print mmCount,line,
