# get correlations between mod frequencies of different cell types
# for the Wang/Xu and Doench experiments
from scipy.stats import linregress, pearsonr, spearmanr, mannwhitneyu, rankdata
from annotateOffs import *
from os.path import basename
import itertools
import glob

def parseTab(fname):
    " return modFreqs of a tab file as a dict seq -> freq "
    print "parsing %s" % fname
    freqs = {}
    for row in iterTsvRows(open(fname)):
        seq = row.seq
        freqs[seq] = float(row.modFreq), float(row.fusi), float(row.wang)
    return freqs

def toLists(d1, d2):
    " give to dictionaries with seq -> freq, convert to two lists with frequencies "
    l1, l2 = [], []
    l1b, l2b = [], []
    l1c, l2c = [], []
    allSeqs = set(d1.keys()).intersection(d2.keys())
    for seq in allSeqs:
        l1.append(d1[seq][0])
        l2.append(d2[seq][0])
        l1b.append(d1[seq][1])
        l2b.append(d2[seq][1])
        l1c.append(d1[seq][2])
        l2c.append(d2[seq][2])
    return l1, l2, l1b, l2b, l1c, l2c

def compareFreqs(fname1, fname2, ofh, neg=False):
    f1 = parseTab(fname1.replace(".scores.tab","")+".scores.tab")
    f2 = parseTab(fname2.replace(".scores.tab","")+".scores.tab")
    l1, l2, l1fusi, l2fusi, l1wang, l2wang = toLists(f1, f2)
    if len(l1)==0:
        return
    if neg:
        l1 = [-x for x in l1]
        l2 = [-x for x in l2]
    print "comparing %s, %s, lengths: %d, %d" % (fname1, fname2, len(l1), len(l2))
    #pearR, pearP = pearsonr(l1, l2)
    #fusiR1, fusiP1 = pearsonr(l1, l1fusi)
    #fusiR2, fusiP2 = pearsonr(l2, l2fusi)
    #wangR1, wangP1 = pearsonr(l1, l1wang)
    #wangR2, wangP2 = pearsonr(l2, l2wang)
    pearR, pearP = spearmanr(l1, l2)
    fusiR1, fusiP1 = spearmanr(l1, l1fusi)
    fusiR2, fusiP2 = spearmanr(l2, l2fusi)
    wangR1, wangP1 = spearmanr(l1, l1wang)
    wangR2, wangP2 = spearmanr(l2, l2wang)
    #spearR, spearP = spearmanr(l1, l2)
    #print "%s %s: pearson R %0.3f (p %0.3f)" % (basename(fname1), basename(fname2), pearR, pearP)

    fname1= fname1.replace("Train", "-")
    fname2= fname2.replace("Train", "-")
    study = basename(fname1).split(".")[0].split("-")[0]
    cell1 = basename(fname1).split(".")[0].split("-")[-1]
    cell2 = basename(fname2).split(".")[0].split("-")[-1]
    #row = [study, cell1, cell2, "%0.3f" % pearR, "%0.3f" % spearR]
    row = [study, cell1, cell2, "%0.3f" % pearR, "%0.3f" % fusiR1, "%0.3f" % fusiR2, "%0.3f" % wangR1, "%0.3f" % wangR2]
    ofh.write("\t".join(row)+"\n")

ofh = open("out/corrCellType.tsv", "w")
#headers = ["study", "cellType1", "cellType2", "pearsonR", "spearmanR"]
headers = ["study", "cellType1", "cellType2", "modFreqSpearmanR", "fusiModFreq1SpearmanR", "fusiModFreq2SpearmanR", "wangModFreqSpearR1", "wangModFreqSpearR2"]
ofh.write("\t".join(headers)+"\n")

compareFreqs("effData/xu2015TrainHl60", "effData/xu2015TrainKbm7", ofh, neg=True)
compareFreqs("effData/xu2015TrainMEsc1", "effData/xu2015TrainMEsc2", ofh, neg=True)
#compareFreqs("effData/xu2015TrainMEsc1.guides.tab", "effData/xu2015TrainMEsc.guides.tab", ofh)
#compareFreqs("effData/xu2015TrainMEsc2.guides.tab", "effData/xu2015TrainMEsc.guides.tab", ofh)
compareFreqs("effData/doench2014-Hs-TF1_CD13", "effData/doench2014-Hs-NB4_CD13", ofh)
compareFreqs("effData/doench2014-Hs-TF1_CD33", "effData/doench2014-Hs-NB4_CD33", ofh)
compareFreqs("effData/doench2014-Hs-MOLM13_CD33", "effData/doench2014-Hs-TF1_CD33", ofh)
compareFreqs("effData/doench2014-Hs-NB4_CD33", "effData/doench2014-Hs-MOLM13_CD33", ofh)
compareFreqs("effData/chari2015Train293T", "effData/chari2015TrainK562", ofh)

hartFnames = glob.glob("effData/hart2016-*Avg.scores.tab")
#hartFnames = glob.glob("effData/hart2016*.scores.tab")
#hartFnames = [fn for fn in hartFnames if "Pseudo" not in fn and "Avg" in fn]
#hartFnames = [fn for fn in hartFnames if "Avg" in fn]
for fn1, fn2 in itertools.combinations(hartFnames, 2):
    print fn1, fn2
    compareFreqs(fn1, fn2, ofh)


ofh.close()
print "wrote output to %s" % (ofh.name)

