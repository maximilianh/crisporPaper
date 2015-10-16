# plot GC content and off-target count as a scatter plot

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
import matplotlib.pyplot as plt
import numpy
import matplotlib.backends.backend_pdf as pltBack
import operator
from matplotlib import gridspec
import string
from annotateOffs import *

from collections import defaultdict

def parseOtRatios():
    otCounts = defaultdict(int)
    gcCounts = defaultdict(int)
    allGuides = set()
    offtargetFreqs = defaultdict(float)
    ontargetFreqs = dict()
    guideSeqs = dict()
    for row in iterTsvRows("offtargets.tsv"):
        if row.type==("on-target"):
            guideGc = (row.seq[:20].count("G") + row.seq[:20].count("C")) / 20.0
            if row.score=="NA":
                continue
            gcCounts[row.name] = guideGc
            ontargetFreqs[row.name] = float(row.score)
            guideSeqs[row.name] = row.seq
            continue
        allGuides.add(row.name)
        freq = float(row.score)
        if freq==0.0:
            continue
        #if freq<0.01:
            #print "skipping", row
            #continue
        otCounts[row.name]+=1
        offtargetFreqs[row.name] += freq

    offtargetFreqs["Ran_EMX1-sg1"] = 0.0
    offtargetFreqs["Ran_EMX1-sg2"] = 0.0

    #for row in iterTsvRows("origData/Ran2015/convertNoSeq.tab"):
        #freq = float(row.score)
        #if row.type==("on-target"):
            #ontargetFreqs[row.name] = freq
            #continue
        #print row.name, freq
        #offtargetFreqs[row.name] += freq

    #print "ontarget", ontargetFreqs["Ran_EMX1-sg2"]
    #print "offtarget", offtargetFreqs["Ran_EMX1-sg2"]

    offtargetRatios = {}
    for name, ontargetFreq in ontargetFreqs.iteritems():
        offtargetFreq = offtargetFreqs[name]
        ratio = offtargetFreq / ontargetFreq
        #ratio = ontargetFreq / offtargetFreq
        #if offtargetFreq==0.0:
            #continue

        offtargetRatios[name] = ratio
    print "missing from plot: %s" % (allGuides - set(offtargetRatios))
    print "total number of guides used: %d" % len(offtargetRatios)
    print
    return gcCounts, offtargetRatios

gcCounts, otRatios = parseOtRatios()
#plotData = offtargetRatios

rows = gcCounts.items()
#for name, gcCount in sorted(rows, key=operator.itemgetter(1)):
    #print name, gcCount, guideSeqs[name]

studyX = defaultdict(list)
studyY = defaultdict(list)
#studyZ = defaultdict(list)
studyGuides = defaultdict(list)
for name in otRatios:
    gcCount = gcCounts[name]
    study = name.split("_")[0]
    yVal = otRatios[name]
    studyX[study].append(gcCount)
    studyY[study].append(yVal)
    studyGuides[study].append(name)

colors  = ["green", "blue", "green", "blue", "red", "grey", "orange", "blue"]
markers = ["o", "s", "+", ">", "<", "o", ".", "o"]
figs = []
i = 0
studyNames = []

gs  = gridspec.GridSpec(2, 1, height_ratios=[1, 5])
ax = plt.subplot(gs[0])
ax2 = plt.subplot(gs[1])
#f,(ax,ax2) = plt.subplots(2,1,sharex=True)

studies = studyX.keys()
studies.sort()

for study in studies:
    xVals = studyX[study]
    yVals = studyY[study]
    guideNames = studyGuides[study]
    #zVals = studyZ[study]
    for a in [ax, ax2]:
        studyFig = a.scatter(xVals, yVals, \
            alpha=0.9, \
            marker=markers[i], \
            s=30, \
            color=colors[i])
    figs.append(studyFig)
    studyNames.append(study)
    i+=1

    for x, y, guideName in zip(xVals, yVals, guideNames):
           # arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0')
           # bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5))
           guideName = string.split(guideName, "_", maxsplit=1)[1]
           if float(x)>=0.75 or y > 2:
               for a in [ax, ax2]:
                   a.annotate(
                      guideName, fontsize=9, rotation=0, ha="right", rotation_mode="anchor",
                      xy = (x, y), xytext = (-7,0), alpha=1.0,
                      textcoords = 'offset points', va = 'bottom')

plt.legend(figs,
       studyNames,
       scatterpoints=1,
       loc='upper left',
       ncol=1,
       fontsize=10)

ax.set_ylim(20,55) # outliers only
ax2.set_ylim(-.5,10) # most of the data

# from http://matplotlib.org/examples/pylab_examples/broken_axis.html
ax.spines['bottom'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax.xaxis.tick_top()
ax.tick_params(labeltop='off') # don't put tick labels at the top
ax2.xaxis.tick_bottom()
ax.set_yticks([20,30,40,50])

d = .015 # how big to make the diagonal lines in axes coordinates
# arguments to pass plot, just so we don't keep repeating them
kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
ax.plot((-d,+d),(-d,+d), **kwargs)      # top-left diagonal
ax.plot((1-d,1+d),(-d,+d), **kwargs)    # top-right diagonal

kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
ax2.plot((-d,+d),(1-d,1+d), **kwargs)   # bottom-left diagonal
ax2.plot((1-d,1+d),(1-d,1+d), **kwargs) # bottom-right diagonal

plt.subplots_adjust(hspace=0.15)

plt.xlabel("GC content of guide sequence")
#plt.ylabel("Number of off-targets with mod. freq. > 0.1%")
#plt.ylabel("Number of off-targets")
plt.ylabel("Ratio total off-target / on-target frequency")
#plt.ylim(0,50)
outfname = "out/gcOtCount.pdf"
plt.savefig(outfname, format = 'pdf')
plt.savefig(outfname.replace(".pdf", ".png"))
print "wrote out/gcOtCount.pdf / .png"

