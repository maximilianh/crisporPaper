# run all guide sequences through both Fusi implementations and compare both results

from annotateOffs import *
from scipy.stats import pearsonr, spearmanr, mannwhitneyu, rankdata
sys.path.insert(0, "../crispor/")
from crisporEffScores import *

setBinDir("../crispor/bin")

def compareFusiImplementations():
    """
    run all guide sequences through both Fusi implementations and compare them
    """
    oldScores = []
    seqs = []
    for fname in glob.glob("effData/*.scores.tab"):
        print "parsing", fname
        for row in iterTsvRows(fname):
            oldScore = int(round(float(row.fusi)*100))
            oldScores.append(oldScore)
            seq = trimSeqs([row.longSeq100Bp], -24, 6)
            seqs.append(seq[0])
            #newScore = calcFusiDoench()

    seqToOldScore = {}
    for seq, score in zip(seqs, oldScores):
        seqToOldScore[seq] = score
        
    print "running fusi"
    uniqueSeqs = list(set(seqs))
    newScores = calcFusiDoench(uniqueSeqs)

    seqToNewScore = {}
    for seq, score in zip(uniqueSeqs, newScores):
        seqToNewScore[seq] = score

    ofh = open("out/fusiOldNew.tab", "w")
    news = []
    olds = []
    for seq, oldScore in seqToOldScore.iteritems():
        newScore = seqToNewScore[seq]
        #print oldScore, newScore
        ofh.write("%s\t%d\t%d\n" % (seq, oldScore, newScore))
        news.append(newScore)
        olds.append(oldScore)

    print "spearman", pearsonr(olds, news)
    print "wrote values to %s" % ofh.name

compareFusiImplementations()
