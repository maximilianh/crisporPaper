import random, string

# plot on-target read fraction versus Gc content
# HAS A PARAMETER THAT CAN PLOT RELATIVE FREQUENCIES!!

from annotateOffs import *

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.backends.backend_pdf as pltBack

#ignoreStudies = ["Frock", "Cho"]
#ignoreStudies = ["Cho"]
ignoreStudies = []

plotRel = False

def plot(maxMMs, guideGcs, otCounts):
    " maxMMs is a dict guide name -> maximum mismatch  "
    xVals = []
    yVals = []
    zVals = []

    studyX = defaultdict(list)
    studyY = defaultdict(list)
    studyZ = defaultdict(list)
    studyLabels = defaultdict(list)
    for name, maxMM in maxMMs.iteritems():
        study = name.split("_")[0]
        xVals.append(guideGcs[name])
        studyX[study].append( guideGcs[name] )
        studyY[study].append( maxMM )
        studyZ[study].append( otCounts[name] )
        guideName = string.split(name, "_", 1)[1]

        studyLabels[study].append(guideName )

    outfname = "gcReadFraction" + '.pdf'
    pdf = pltBack.PdfPages(outfname)
    fig = plt.figure(figsize=(5,5),
                   dpi=300, facecolor='w')
    fig = plt.figure()

    colors = ["red", "blue", "black", "black", "green", "grey", "orange", "violet"]
    markers = ["o", "s", "<", "+", "+", "^", ".", ">"]
    studyNames = []
    figs = []
    i = 0
    studies = sorted(studyX.keys())
    print studies
    for study in studies:
        xVals = studyX[study]
        yVals = studyY[study]
        zVals = studyZ[study]
        #plt.scatter(xVals, yVals, \
            #alpha=.8, \
            #marker="o", \
            #s=zVals, \
            #color="grey")
        studyFig = plt.scatter(xVals, yVals, \
            alpha=.5, \
            marker=markers[i], \
            s=30, \
            color=colors[i])
        figs.append(studyFig)
        studyNames.append(study)
        i+=1

        labels = studyLabels[study]
        for x, y, label in zip(xVals, yVals, labels):
               # arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0')
               # bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5))
               plt.annotate(
                  label, fontsize=8, rotation=30, ha="right", rotation_mode="anchor",
                  xy = (x, y), xytext = (0,0), alpha=0.9,
                  textcoords = 'offset points', va = 'bottom')

    plt.legend(figs,
           studyNames,
           scatterpoints=1,
           loc='lower left',
           ncol=2,
           fontsize=10)

    plt.ylim(0,1.08)
    plt.xlabel("GC content")
    if plotRel:
        plt.ylabel("relative efficacy (fraction of reads relative to all indel-causing reads)")
    else:
        plt.ylabel("On-target modification frequency")
    fig.savefig(pdf, format = 'pdf')
    print "Wrote %s" % outfname

    outfname = outfname.replace(".pdf", ".png")
    fig.savefig(outfname)

    pdf.close()
    print "Wrote %s" % outfname

def annotateOts():
    " write annotations of the offtargets to a tab-setp file and also return as dict "
    targetSeqs = {}
    inFname = "offtargets.tsv"
    inRows = []
    otCounts = defaultdict(int)
    for row in iterTsvRows(inFname):
        study = row.name.split("_")[0]
        if study in ignoreStudies:
            continue
        if row.type=="on-target":
            targetSeqs[row.name] = row.seq
        else:
            otCounts[row.name] += 1
        inRows.append(row)

    # first sum up the frequencies for each guide
    sums = defaultdict(float)
    for row in inRows:
        if row.score=="NA": # Frock cannot quantify the target
            continue
        sums[row.name] += float(row.score)

    headers = ["name", "guideSeq", "otSeq", "guideGc", "otGc", "assayScore", "mmCount", "diffLogo", "mmCountOneGap", "oneGapSeq", "diffLogoOneGap"]
    rows = [headers]

    maxMMs = defaultdict(int)
    guideGcConts = {}
    i = 0
    for row in inRows:
        #if row.type=="on-target":
            #continue
        guideSeq = targetSeqs[row.name][:-3].upper()
        otSeq = row.seq[:-3].upper()
        mmCount, diffLogo = countMmsAndLogo(guideSeq, otSeq)

        otGcCont = gcCont(otSeq)
        guideGc = gcCont(guideSeq[:20])

        gappedMm, guideGapSeqs, otGapSeqs, gapLogos = findGappedSeqs(guideSeq, otSeq)
        if gappedMm > mmCount:
            gappedMm = 0
            guideGapSeqs = []
            gapLogos = []

        otRow = [row.name, guideSeq, otSeq, str(guideGc), str(otGcCont), row.score, str(mmCount), diffLogo, str(gappedMm), ",".join(guideGapSeqs), ",".join(gapLogos)]
        rows.append(otRow)
        #dataName = row.name+str(i)
        dataName = row.name
        #maxMMs[row.name+str(i)] = max(maxMMs[row.name], mmCount)
        #maxMMs[dataName] = mmCount + random.random()
        #maxMMs[dataName] = mmCount
        if row.score=="NA": # Frock cannot quantify the target
            continue
        freq = float(row.score)
        if plotRel:
            if "Tsai" not in dataName:
                freq = freq/sums[dataName]
                #print "correction: ", dataName, freq, sums[dataName]
        #else:
            #freq = 1.0 - freq
        if row.type=="off-target":
            continue
        assert(dataName not in maxMMs)
        maxMMs[dataName] = freq

        #guideGcConts[row.name] = guideGc
        #guideGcConts[dataName] = otGcCont+random.random()*3
        guideGcConts[dataName] = guideGc
        i+=1

    ofh = open("annotOfftargets.tsv", "w")
    for row in rows:
        ofh.write("\t".join(row))
        ofh.write("\n")
    ofh.close()
    print "wrote %s" % ofh.name

    return maxMMs, guideGcConts, otCounts

def main():
    maxMMs, guideGcConts, otCounts = annotateOts()
    plot(maxMMs, guideGcConts, otCounts)

main()
