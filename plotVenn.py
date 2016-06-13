# plot the overlaps for the guides VEGFA and EMX1 for which we have multiple studies

# note to self: use homebrew's python to run this /usr/local/bin/python2.7 

import os, logging
from annotateOffs import *
from collections import defaultdict
#from scipy.stats import linregress

import matplotlib as mpl
mpl.use('Agg')

from matplotlib_venn import venn3, venn2

import matplotlib.pyplot as plt
import numpy as np

def readMitPredOts(fnames):
    " return a set with all MIT predicted off-targets "
    ret = set()
    for fname in fnames:
        for line in open(fname):
            fs = line.split(", ")
            ret.add(fs[4])
    return ret

def parseAllOffs(fname, skipCell=False):
    """ return a dict targetSeq -> study -> list of (off-target sequence, modFreq)
    and a dict targetSeq -> targetName
    """
    nameToTarget = defaultdict(dict)
    targetToName = {}
    for row in iterTsvRows(fname):
        if row.type=="on-target":
            nameToTarget[row.name] = row.seq
            targetToName[row.seq] = row.name

    seqs = defaultdict(dict)
    for row in iterTsvRows(fname):
        if row.type!="off-target":
            continue
        targetSeq = nameToTarget[row.name]
        study = row.name.split("_")[0]
        if skipCell:
            study = study.split("/")[0]
        modFreq = float(row.score)
        #if modFreq<0.001:
            #continue
        seqs[targetSeq].setdefault(study, []).append( (row.seq, modFreq) ) 

    # special filtering: make two indexes into data to facilitate removal
    studyCounts = defaultdict(int) # offtarget -> number of studies
    maxFreq = defaultdict(float) # offtarget -> highest mod freq
    for targetSeq, studyFreqs in seqs.iteritems():
        for study, seqFreqs in studyFreqs.iteritems():
            for seq, freq in seqFreqs:
                studyCounts[seq] +=1
                maxFreq[seq] = max(maxFreq[seq], freq)

    # special filtering: remove offtarget < 0.001 when found only by one study
    # and with a freq < 0.001
    remOts = set()
    for otSeq, count in studyCounts.iteritems():
        if count != 1:
            continue
        if maxFreq[otSeq] < 0.001:
            remOts.add(otSeq)

    filtSeqs = defaultdict(dict)
    for targetSeq, studyFreqs in seqs.iteritems():
        for study, seqFreqs in studyFreqs.iteritems():
            for seq, freq in seqFreqs:
                if seq in remOts:
                    continue
                filtSeqs[targetSeq].setdefault(study, []).append((seq, freq))

    #assert("Kim16" in filtSeqs["GGGTGGGGGGAGTTTGCTCCTGG"])
    return filtSeqs, targetToName

def createTsv(overlapGuides):
    # create the TSV files
    outTsvFnames = []
    seqs, seqNames = parseAllOffs("offtargets.tsv")
    mitSeqs = readMitPredOts(["mitOfftargets/Hsu_EMX1.3.csv", "mitOfftargets/Frock_VEGFA.csv"])

    for guideSeq, guideName in overlapGuides:
        # prep for tsv file: convert to dict otSeq -> study -> freq
        studySeqs = seqs[guideSeq]

        otFreqs = defaultdict(dict)
        studies = sorted(studySeqs)
        #print guideName, studies, studySeqs
        for studyName in studies:
            otInfo = studySeqs[studyName]
            for otSeq, otFreq in otInfo:
                otFreqs[otSeq][studyName] = otFreq

        # output to tsv file
        ofh = open("out/venn-%s.tsv" % guideName, "w")
        headers = ["guide", "off-target", "mismatch", "diffLocs", "offtScore"]
        headers.extend(studies)
        headers.append("mitPredicted")
        ofh.write("\t".join(headers)+"\n")

        rows = []
        for otSeq, studyFreqs in otFreqs.iteritems():
            mmCount, diffLogo = countMms(guideSeq, otSeq)
            otScore = calcHitScore(guideSeq, otSeq)
            isMitPred = (otSeq in mitSeqs)
            row = [guideSeq, otSeq, mmCount, diffLogo, otScore]
            for study in studies:
                freq = studyFreqs.get(study, "notFound")
                if study=="Hsu" and freq=="notFound":
                    freq = "notTested"
                row.append(str(freq))
            row.append(str(isMitPred))
            rows.append(row)

        rows.sort(key=operator.itemgetter(2))
        for row in rows:
            ofh.write("\t".join([str(x) for x in row])+"\n")

        outTsvFnames.append(ofh.name)

    print "wrote data to %s" % ", ".join(outTsvFnames)

def main():
    #plt.figure(figsize=(8,10))

    seqs, seqNames = parseAllOffs("offtargets.tsv", skipCell=True)
    
    #overlapGuides = [
        #("GAGTCCGAGCAGAAGAAGAAGGG", "EMX1"),
        #("GGGTGGGGGGAGTTTGCTCCTGG", "VEGFA")
        #]
    overlapGuides = []
    for targetSeq, studyOffs in seqs.iteritems():
        #print targetSeq, studyOffs.keys(), seqNames[targetSeq]
        parts = seqNames[targetSeq].split("_")
        if len(parts)==2:
            name = parts[-1]
        if len(parts)==3:
            name = parts[1]+" "+parts[2]
        if len(studyOffs)>=2:
            overlapGuides.append( (targetSeq, name) )
    #print overlapGuides

    #plotCount = plt.subplots(len(overlapGuides), 1)
    plotCount = 2
    fig, axArr = plt.subplots(plotCount, 1)
    fig.set_size_inches(5,plotCount*5)

    studyDescs = {
        "Tsai" : "GuideSeq\n(Tsai et al, HEK293)",
        "Frock" : "Translocation\nsequencing\n(HTGTS,\nFrock et al, HEK293T)",
        "Hsu" : "Targeted sequencing\n(Hsu et al,\nHEK293FT)",
        "Kim" : "DigenomeSeq\n(Kim et al 2015, HAP1)",
        "Kim16" : "DigenomeSeq2\n(Kim et al 2016, HeLa)",
    }

    plotRow = 0
    for (guideSeq, guideName) in overlapGuides:
        if guideName not in ["EMX1", "VEGFA site1"]:
            continue

        print "guide: ",guideName
        studySeqs = seqs[guideSeq]

        ax = axArr[plotRow]
        labels = []
        sets = []
        for studyName, seqInfo in studySeqs.items():
            #if studyName not in studyDescs:
                #continue
            if studyName=="Kim" and guideName=="EMX1":
                continue
            if studyName=="Hsu":
                continue
            print "using guides from ", studyName
            studyName = studyDescs[studyName]
            labels.append(studyName)
            studyOts = set()
            for otSeq, otFreq in seqInfo:
                if float(otFreq)!=0.0:
                    studyOts.add(otSeq)
            print "offtarget-count:", len(studyOts)
            sets.append(set(studyOts))

        if len(sets)==3:
            #print guideName, sets, labels, ax
            venn3(subsets=sets, set_labels=labels, ax=ax, labelSize="small")
            ax.set_title(guideName)
        elif len(sets)==2:
            venn2(subsets=sets, set_labels=labels, ax=ax)
            ax.set_title(guideName)
        else:
            print len(sets)
            assert(False)

        plotRow += 1


    fig.tight_layout()

    outFname = "out/venn.pdf"
    plt.savefig(outFname)
    plt.savefig(outFname.replace(".pdf", ".png"))
    print "wrote plot to %s and .png" % outFname

    createTsv(overlapGuides)

main()
