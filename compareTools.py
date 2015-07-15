# compare MIT crispr site with other sites

# first parameter: the off-target score cutoff for crispor

import annotateOffs
import glob, sys
from os.path import basename, splitext
from collections import defaultdict

#ignoreStudies = ["Hsu", "Cho"]
#ignoreStudies = ["Hsu", "Cho"]
ignoreStudies = []

def parseGoldStandard(fname, minFrac):
    " return all curated off-targets with a minimum score of minFrac "
    targets = dict()
    data = defaultdict(set)
    for line in open(fname).read().splitlines()[1:]:
        fs = line.split("\t")
        sType = fs[-1]
        guideName = fs[0]
        study = guideName.split("_")[0]
        if study in ignoreStudies:
            continue
        seq = fs[1]
        if sType=="on-target":
            targets[guideName] = seq
            continue
        frac = fs[2]
        if float(frac) < minFrac:
            continue
        data[guideName].add(seq)
    return targets, data


def parseCrispor(dirName, minOtScore, minAltOtScore, targets):
    " returns a tuple of (dict guideName -> on-target seq, dict guideName -> set of offtargetSeqs) "
    fnames= glob.glob(dirName+"/*.tsv")
    data = defaultdict(set)
    for fname in fnames:
        guideName = splitext(basename(fname))[0]
        study = guideName.split("_")[0]
        if study in ignoreStudies:
            continue
        if guideName not in targets:
            print ("no benchmark data for crispor data for %s" % guideName)
            continue
        for line in open(fname):
            if line.startswith("guideId"):
                continue
            fs = line.strip().split("\t")
            targetSeq = fs[1]
            otSeq = fs[2]
            otScore = fs[4]
            if otSeq.endswith("GA") or otSeq.endswith("AG"):
                if float(otScore) < minAltOtScore:
                    continue
            else:
                if float(otScore) < minOtScore:
                    continue

            data[guideName].add(otSeq)
    return data

def parseMit(dirName, targetSeqs):
    " parse the MIT csv files, return a dict with base-filename -> set of offtargets "
    fnames= glob.glob(dirName+"/*.csv")
    data = defaultdict(set)
    for fname in fnames:
        guideName = splitext(basename(fname))[0]
        study = guideName.split("_")[0]
        if study in ignoreStudies:
            continue
        if guideName not in targetSeqs:
            print "no benchmark data for MIT off-target data %s" % guideName
            continue
        for line in open(fname):
            if line.startswith("guide"):
                continue
            fs = line.split(", ")
            seq = fs[4]
            if seq == targetSeqs[guideName]:
                continue
            data[guideName].add(seq)
    return data

def compareOffs(benchSeqs, names, offList):
    ofh = open("out/toolComparison.tsv", "w")
    offData = {}
    headers = ["guide", "validOffs"]
    for name in names:
        headers.append("Count_%s" % name)
        headers.append("Prec_%s" % name)
        headers.append("Recall_%s" % name)
    ofh.write("\t".join(headers)+"\n")

    # global set across all studies of all predictions and all validated off-targets
    allPreds = defaultdict(set) # dict tool -> predicted seqs
    allBench = set() # set of all validated off-targets

    for guideName, benchOffSeqs in benchSeqs.iteritems():
        allBench.update(benchOffSeqs)
        row = [guideName, str(len(benchOffSeqs))]

        for toolName, offDict in zip(names, offList):
            offSeqs = offDict[guideName]
            allPreds[toolName].update(offSeqs)

            correctPred = len(benchOffSeqs.intersection(offSeqs))

            #print "corrPred, offSeqs", correctPred, len(offSeqs)
            if len(offSeqs)!=0:
                prec = float(correctPred) / len(offSeqs)
            else:
                prec = 0.0

            if len(benchOffSeqs)!=0:
                recall = correctPred / float(len(benchOffSeqs))
            else:
                recall = 0.0

            row.append ("%d" % len(offSeqs))
            row.append ("%0.3f" % prec)
            row.append ("%0.2f" % recall)
        ofh.write("\t".join(row)+"\n")
    ofh.close()
    print "Wrote per-guide results to %s" % ofh.name

    print "Total off-targets to find: %d" % len(allBench)
    headers = ["tool", "totalPred", "prec", "recall"]
    for toolName, predSet in allPreds.iteritems():
        correctPreds = predSet.intersection(allBench)
        prec = float(len(correctPreds)) / len(predSet)
        recall = float(len(correctPreds)) / len(allBench)
        row = [toolName, str(len(predSet)), "%0.2f" % prec, "%0.2f" % recall]
        print "\t".join(row)

if len(sys.argv)==1:
    print "WARNING: no arguments specified, using default values minFrac=0.01, minOtScore=0.1 and minAltOtScore=0.2"
    minFrac = 0.01
    minOtScore = 0.1
    minAltOtScore = 0.1
else:
    minFrac = float(sys.argv[1])
    minOtScore = float(sys.argv[2])
    minAltOtScore = float(sys.argv[3])

targets, benchOts = parseGoldStandard("offtargetsFilt.tsv", minFrac)

crisporOffs = parseCrispor("crisporOfftargets/", minOtScore, minAltOtScore, targets)
mitOffs = parseMit("mitOfftargets/", targets)

# make sure we have OT data for all guides
benchNotMit = set(benchOts) - set(mitOffs)
if len(benchNotMit)!=0:
    print "Guides missing in MIT: ", benchNotMit
benchNotCrispor = set(benchOts) - set(crisporOffs)
if len(benchNotCrispor)!=0:
    print "Guides Missing in Crispor: ", benchNotCrispor


compareOffs(benchOts, ["CRISPOR", "MIT"], [crisporOffs, mitOffs])

