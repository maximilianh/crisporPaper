# run all guide sequences through both wang implementations and compare both results

from annotateOffs import *
from crisporEffScores import *
from scipy.stats import pearsonr, spearmanr, mannwhitneyu, rankdata
import crisporEffScores

setBinDir("../crispor/bin")

def compareWangImplementations():
    """
    run all guide sequences through both wang implementations and compare them
    """
    print "parsing"
    seqs, scores = parseAllGuides()
    seqs = [s[:20].upper() for s in seqs]
    seqs = list(set(seqs))
    print "using libsvm directly"
    newScores = calcWangSvmScores(seqs)
    print "using R"
    origScoreDict = calcWangSvmScoresUsingR(seqs)

    origScores = []
    for s in seqs:
        origScores.append(origScoreDict[s])

    ofh = open("out/wangDiffs.tsv", "w")
    ofh.write("seq\torig\tnew\tdiff\n")
    for seq, origScore, newScore in zip(seqs, origScores, newScores):
        row = [seq, origScore, newScore, newScore-origScore]
        row = [str(x) for x in row]
        ofh.write("\t".join(row))
        ofh.write("\n")

    print "Pearson correlation and p-Value", pearsonr(origScores, newScores)
    print "details written to %s" % ofh.name

compareWangImplementations()
