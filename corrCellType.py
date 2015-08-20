# get correlations between mod frequencies of different cell types
# for the Xu and Doench experiments
from scipy.stats import linregress, pearsonr, spearmanr, mannwhitneyu, rankdata
from annotateOffs import *
from os.path import basename

def parseTab(fname):
    " return modFreqs of a tab file as a dict seq -> freq "
    freqs = {}
    for row in iterTsvRows(open(fname)):
        seq = row.seq
        freqs[seq] = float(row.modFreq)
    return freqs

def toLists(d1, d2):
    " give to dictionaries with seq -> freq, convert to two lists with frequencies "
    #assert(len(d1)==len(d2))
    l1, l2 = [], []
    for seq in d1:
        l1.append(d1[seq])
        l2.append(d2[seq])
    #print len(d1), len(d2)
    return l1, l2

def compareFreqs(fname1, fname2, ofh):
    f1 = parseTab(fname1)
    f2 = parseTab(fname2)
    l1, l2 = toLists(f1, f2)
    pearR, pearP = pearsonr(l1, l2)
    #print "%s %s: pearson R %0.3f (p %0.3f)" % (basename(fname1), basename(fname2), pearR, pearP)
    fname1= fname1.replace("Train", "-")
    fname2= fname2.replace("Train", "-")
    study = basename(fname1).split(".")[0].split("-")[0]
    cell1 = basename(fname1).split(".")[0].split("-")[-1]
    cell2 = basename(fname2).split(".")[0].split("-")[-1]
    row = [study, cell1, cell2, "%0.3f" % pearR]
    ofh.write("\t".join(row)+"\n")

ofh = open("out/corrCellType.tsv", "w")
headers = ["study", "cellType1", "cellType2", "pearsonR"]
ofh.write("\t".join(headers)+"\n")

compareFreqs("effData/xu2015TrainHl60.ext.tab", "effData/xu2015TrainKbm7.ext.tab", ofh)
compareFreqs("effData/doench2014-Hs-TF1_CD13.ext.tab", "effData/doench2014-Hs-NB4_CD13.ext.tab", ofh)
compareFreqs("effData/doench2014-Hs-TF1_CD33.ext.tab", "effData/doench2014-Hs-NB4_CD33.ext.tab", ofh)
compareFreqs("effData/doench2014-Hs-MOLM13_CD33.ext.tab", "effData/doench2014-Hs-TF1_CD33.ext.tab", ofh)
compareFreqs("effData/doench2014-Hs-NB4_CD33.ext.tab", "effData/doench2014-Hs-MOLM13_CD33.ext.tab", ofh)
compareFreqs("effData/chari2015Train293T.tab", "effData/chari2015TrainK562.tab", ofh)

ofh.close()
print "wrote output to %s" % (ofh.name)

