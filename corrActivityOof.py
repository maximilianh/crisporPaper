# determine correlation between guide activity and the OOF score

from annotateOffs import *
import os
import glob
from os.path import basename, splitext
from scipy.stats import linregress, pearsonr, spearmanr
import operator

import matplotlib
matplotlib.use('Agg')
#matplotlib.rcParams['pdf.fonttype'] = 42
import matplotlib.pyplot as plt
import numpy as np
from numpy import mean

rows = []

fig, axArr = plt.subplots( 40, 1 )
fig.set_size_inches(5,30)

axIdx = 0
for fname in glob.glob("effData/*.tab"):
    if ".ext." in fname:
        continue
    if "museum" in fname:
        continue
    if "doench2014-Mm" in fname:
        print "XX ignoring Mouse data for now, too many duplicates"
        continue
    if "CD33Exon2" in fname or "xu2015." in fname:
        # doench-single exon is too small
        continue
    if "Valid" in fname:
        # skip small-scale data, not relevant here
        continue
    
    oofScores = []
    actScores = []
    for row in iterTsvRows(fname):
        seq = row.seq
        if seq=="GGGAGGAGATAAGAAGAGAAAGG":
            # one sequence in Doench has a 1bp deletion relative to hg19
            continue
        if seq=="TGGACATCTTGGGCGTGGTGTGG":
            # one sequence wasn't mapped in Doench
            continue
        if seq=="GCAAGGATCTTCAAAAAGCACGG":
            # one farboud sequence wasn't mapped
            continue
        if seq=="GGGACTCACTGGGAACGGTGTGG" or seq=="GGGACTCACTGGGAACGGTGTGG" or seq=="GGGACTCACTGGGAACGGTGTGG" or seq=="GGGACTCACTGGGAACGGTGTGG" or seq=="GGAGATCGGCTCCCTTCTATTGG" or seq=="AGACCATCGCCGCTTCACTAAGG" or seq=="GACATTTGCAGAAGATGAAGTGG" or seq=="AGTCGTGAGTTTGGAGCGACGGG":
            # one gagnon sequence
            continue
        oofScore = lookupOofScore(seq)
        if oofScore is None:
            print "no oof score for %s" % seq
            continue
        oofScores.append(oofScore)
        koScore = row.modFreq
        if "doench" in fname and not "-Hs" in fname:
            koScore = row.geneAbund
        actScores.append(float(koScore))
    pearsonR, pearsonP = pearsonr(actScores, oofScores)
    spearmanR, spearmanP = spearmanr(actScores, oofScores)
    #pearsonR, pearsonP = pearsonr(actScores, oofScores)
    #speTGGGTGGAGGAGTAAGGAATGGGarmanR, spearmanP = spearmanr(actScores, oofScores)
    dataset = basename(fname).split(".")[0]
    rows.append((dataset, pearsonR, pearsonP, spearmanR, spearmanP))

    ax = axArr[axIdx]
    axIdx += 1
    ax.scatter(oofScores, actScores, alpha=0.5, s=5)
    ax.set_title(dataset)
    ax.set_xlabel("OOF score")
    ax.set_xlabel("KO score")
    #print "wrote scatter to %s" % pltFname
    
    data = np.array( (oofScores, actScores), dtype=float )
    #data = np.vstack( (np.array(oofScores), np.array(actScores)))
    #print data.shape
    if actScores==[]:
        continue
    highMin = np.percentile(actScores, 80)
    lowMax = np.percentile(actScores, 20)
    highs = []
    lows = []
    for oof, act in zip(oofScores, actScores):
        if act >= highMin:
            highs.append(oof)
        if act <= lowMax:
            lows.append(oof)
    print "%s - mean OOF top 20%%: %0.3f, mean OOF bottom 20%%: %0.3f" % (dataset, mean(highs), mean(lows))


pltFname = "out/corrActivityOof.pdf"
fig.tight_layout()
plt.savefig(pltFname)
print "Wrote plot to %s" % pltFname

rows.sort(key=operator.itemgetter(1))

for row in rows:
    dataset, pearsonR, pearsonP, spearmanR, spearmanP = row
    if dataset.startswith("schoenig"):
        continue
    if pearsonP > 0.05:
        continue
    print "%s: pearson %0.3f (P %0.3f), spearman %0.3f (P %0.3f)" % (dataset, pearsonR, pearsonP, spearmanR, spearmanP)
