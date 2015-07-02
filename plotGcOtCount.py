# plot GC content and off-target count as a scatter plot

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
import matplotlib.pyplot as plt
import numpy
import matplotlib.backends.backend_pdf as pltBack

import annotateOffs

from collections import defaultdict
#data = numpy.genfromtxt("offtargets.tsv", names=True)
otCounts = defaultdict(int)
gcCounts = defaultdict(int)
for line in open("offtargets.tsv"):
    line = line.strip()
    if line.startswith("name"):
        continue
    if line.endswith("on-target"):
        continue
    fs = line.strip().split()
    freq = fs[2]
    if float(freq)<0.01:
        continue
    name = fs[0]
    seq = fs[1]
    seqGc = annotateOffs.gcCont(seq)
    gcCounts[name] = seqGc
    otCounts[name]+=1

print otCounts
print gcCounts

studyX = defaultdict(list)
studyY = defaultdict(list)
#studyZ = defaultdict(list)
studyGuides = defaultdict(list)
for name, gcCount in gcCounts.iteritems():
    study = name.split("_")[0]
    otCount = otCounts[name]
    studyX[study].append(gcCount)
    studyY[study].append(otCount)
    studyGuides[study].append(name)

colors  = ["green", "blue", "black", "blue", "red", "grey", "orange"]
markers = ["o", "s", "+", ">", "<", "^", "."]
figs = []
i = 0
studyNames = []
for study, xVals in studyX.iteritems():
    yVals = studyY[study]
    guideNames = studyGuides[study]
    #zVals = studyZ[study]
    studyFig = plt.scatter(xVals, yVals, \
        alpha=.5, \
        marker=markers[i], \
        s=30, \
        color=colors[i])
    figs.append(studyFig)
    studyNames.append(study)
    i+=1

    for x, y, guideName in zip(xVals, yVals, guideNames):
           # arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0')
           # bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5))
           if float(x)>65 or y>25:
               plt.annotate(
                  guideName, fontsize=9, rotation=30, ha="right", rotation_mode="anchor",
                  xy = (x, y), xytext = (0,0), alpha=.5,
                  textcoords = 'offset points', va = 'bottom')

plt.legend(figs,
       studyNames,
       scatterpoints=1,
       loc='upper left',
       ncol=2,
       fontsize=10)

plt.xlabel("GC content of guide sequence")
plt.ylabel("Number of off-targets with frequency > 0.001")
outfname = "gcOtCount.pdf"
plt.savefig(outfname, format = 'pdf')
plt.savefig(outfname.replace(".pdf", ".png"))
print "wrote gcOtCount.pdf / .png"

