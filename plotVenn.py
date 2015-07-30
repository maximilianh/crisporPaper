# plot the overlaps for the guides VEGFA and EMX1 for which we have multiple studies

# note to self: use homebrew's python to run this /usr/local/bin/python2.7 

import os, logging
from annotateOffs import *
from collections import defaultdict
#from scipy.stats import linregress

from matplotlib_venn import venn3, venn2

import matplotlib.pyplot as plt
import numpy as np

def parseAllOffs(fname, skipCell=False):
    " return a dict targetSeq -> study -> list of (off-target sequence, modFreq)"
    nameToTarget = defaultdict(dict)
    for row in iterTsvRows(fname):
        if row.type=="on-target":
            nameToTarget[row.name] = row.seq

    seqs = defaultdict(dict)
    for row in iterTsvRows(fname):
        if row.type!="off-target":
            continue
        targetSeq = nameToTarget[row.name]
        study = row.name.split("_")[0]
        if skipCell:
            study = study.split("/")[0]
        seqs[targetSeq].setdefault(study, []).append( (row.seq, float(row.score)) )
    return seqs

def createTsv(overlapGuides):
    # create the TSV files
    outTsvFnames = []
    seqs = parseAllOffs("offtargets.tsv")

    for guideSeq, guideName in overlapGuides:
        # prep for tsv file: convert to dict otSeq -> study -> freq
        studySeqs = seqs[guideSeq]

        otFreqs = defaultdict(dict)
        studies = sorted(studySeqs)
        for studyName in studies:
            otInfo = studySeqs[studyName]
            for otSeq, otFreq in otInfo:
                otFreqs[otSeq][studyName] = otFreq

        # output to tsv file
        ofh = open("out/venn-%s.tsv" % guideName, "w")
        headers = ["guide", "off-target", "mismatch", "diffLocs", "offtScore"]
        headers.extend(studies)
        ofh.write("\t".join(headers)+"\n")

        rows = []
        for otSeq, studyFreqs in otFreqs.iteritems():
            mmCount, diffLogo = countMms(guideSeq, otSeq)
            otScore = calcHitScore(guideSeq, otSeq)
            row = [guideSeq, otSeq, mmCount, diffLogo, otScore]
            for study in studies:
                freq = studyFreqs.get(study, "notFound")
                row.append(str(freq))
            rows.append(row)

        rows.sort(key=operator.itemgetter(2))
        for row in rows:
            ofh.write("\t".join([str(x) for x in row])+"\n")

        outTsvFnames.append(ofh.name)

    print "wrote data to %s" % ", ".join(outTsvFnames)
def main():
    #plt.figure(figsize=(8,10))
    fig, axArr = plt.subplots(2, 1)
    fig.set_size_inches(5,10)

    seqs = parseAllOffs("offtargets.tsv", skipCell=True)
    
    overlapGuides = [
        ("GAGTCCGAGCAGAAGAAGAAGGG", "EMX1"),
        ("GGGTGGGGGGAGTTTGCTCCTGG", "VEGFA")
        ]

    for plotRow, (guideSeq, guideName) in enumerate(overlapGuides):
        studySeqs = seqs[guideSeq]

        ax = axArr[plotRow]
        print plotRow, guideSeq, guideName
        labels = []
        sets = []
        for studyName, seqInfo in studySeqs.items():
            labels.append(studyName)
            studyOts = set()
            for otSeq, otFreq in seqInfo:
                if float(otFreq)!=0.0:
                    studyOts.add(otSeq)
            sets.append(set(studyOts))

        if len(sets)!=3:
            assert(False)

        venn3(subsets=sets, set_labels=labels, ax=ax)
        ax.set_title(guideName)

    fig.tight_layout()

    outFname = "out/venn.pdf"
    plt.savefig(outFname)
    plt.savefig(outFname.replace(".pdf", ".png"))
    print "wrote plot to %s and .png" % outFname

    createTsv(overlapGuides)

main()
