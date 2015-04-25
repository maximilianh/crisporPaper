#0=Chromosome	1=Start	2=End	3=Name	4=GUIDE-Seq Reads	5=Strand	6=Cells	7=Targetsite	8=Target_Sequence	9=Offtarget_Sequence	10=20 bp protospacer # mismatches	11=3 bp PAM # mismatches	12=Mismatch Total
from collections import defaultdict

byGuide = defaultdict(list)

sums = defaultdict(int)

rows = []
for line in open("nbt.3117-S2.txt"):
    if "tru" in line or "Chromosome" in line:
        continue
    fields = line.rstrip("\n").split("\t")
    name = "Tsai_"+fields[7]
    seq = fields[8]
    otSeq = fields[9]
    readCount = fields[4]
    if otSeq[:20]==seq[:20]:
        sType = "on-target"
    else:
        sType = "off-target"
    #byGuide[name].append( (otSeq, readCount, sType) )
    row = [name, otSeq, readCount, sType]
    sums[name] += int(readCount)
    rows.append(row)

for row in rows:
    name, otSeq, readCount, sType = row
    readFreq = float(readCount) / sums[name]
    row = [name, otSeq, str(readFreq), sType]
    print("\t".join(row))
