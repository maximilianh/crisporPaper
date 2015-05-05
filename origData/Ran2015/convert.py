from collections import defaultdict
sums = defaultdict(float)
rows = []
for line in open("suppTable7.txt"):
    fs = line.strip().split("\t")
    name = fs[0]
    if "on-target" in name:
        sType = "on-target"
    else:
        sType = "off-target"
    seq = fs[1]+fs[2][:3]
    freq = fs[5]
    if freq=="N.D.":
        continue
    if "target" in name:
        sType = "on-target"
        guideName = "Ran_"+name.replace("EMX1_target","EMX1-sg")
        name = "target"
    else:
        sType = "off-target"

    sums[guideName] += float(freq)

    row = [guideName, seq, str(float(freq)/100), sType]
    rows.append(row)

for row in rows:
    # no need to normalize, Ran data is already normalized
    #relFreq = str(float(row[2]) / sums[row[0]])
    #newRow = [row[0], row[1], row[2], row[3]]
    print "\t".join(row)
