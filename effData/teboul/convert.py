import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress, pearsonr, spearmanr, mannwhitneyu, rankdata

x = []
y = []
colors = []

ofh = open("../teboulVivo_mm9.guides.tab", "w")
ofh2 = open("../teboulVitro_mm9.guides.tab", "w")

ofh.write("guide\tseq\tmodFreq\n")
ofh2.write("guide\tseq\tmodFreq\n")

pairVivoFreqFh = open("pairVivoFreqs.txt", "w")
singleVivoFreqFh = open("singleVivoFreqs.txt", "w")

guideNo = 0

for line in open("email.tsv"):
    fs = line.rstrip("\n").split("\t")
    if line.startswith("Gene"):
        continue
    # ['Tm6sf2', 'Tm6sf2_#27', 'GTAAATACAGTTCAGAGATG', 'AGG', 'ND', '61', '11', '18.0', '/']
    gene, guide, seq, pam, invitroEff, pupCount, numMutants, invivoEff = fs[:8]
    if len(fs)==9 and fs[8]!='':
        comment = fs[8]
    guide = guide.replace("#", "No")
    if "Ift80" in guide:
        print "found Ift80 guide that was retracted by ordering lab"
        continue

    guide = "guide"+str(guideNo)

    if "pair" in comment:
        pairVivoFreqFh.write(invivoEff+"\n")
    else:
        singleVivoFreqFh.write(invivoEff+"\n")

    fullSeq = seq+pam
    fullSeq = fullSeq.upper()
    if len(fullSeq)!=23:
        print "guide %s is not 23bp long" % guide
        continue

    guideNo += 1

    if "VALUE" not in invivoEff and "pair" not in comment:
        row = [guide, fullSeq, invivoEff]
        ofh.write("\t".join(row))
        ofh.write("\n")

    if "VALUE" not in invitroEff and "ND" not in invitroEff:
        row2 = [guide, fullSeq, invitroEff]
        ofh2.write("\t".join(row2))
        ofh2.write("\n")

    if "ND" not in invitroEff and "VALUE" not in invivoEff:
        x.append(float(invitroEff))
        y.append(float(invivoEff))
        if "pair" in comment:
            colors.append("red")
        else:
            colors.append("blue")

plt.xlabel("in vitro efficiency (blue = single, red = paired)")
plt.ylabel("in vivo efficiency")
plt.scatter(x, y, c=colors)
plt.savefig("scatter.png")
print "saved plot to scatter.png, data for %d guides" % len(x)
pearR, pearP = pearsonr(x, y)
print "Pearson R=%f, P=%f" % (pearR, pearP)

print "wrote %s and %s" % (ofh.name, ofh2.name)
print "wrote %s and %s" % (pairVivoFreqFh.name, singleVivoFreqFh.name)
