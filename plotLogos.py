# plot weblogos of the most efficient and least efficient sequences in the doench,chari and xu/wang
# datasets

from annotateOffs import *
import numpy as np
import os

def parseTab(fname, perc):
    " return all seqs with scores that are higher than the perc percentile "
    xList = []
    yList = []
    seqs = []
    scores = []
    for row in iterTsvRows(open(fname)):
        seq = row.seq[:20]
        seqs.append(seq)
        score = float(row.modFreq)
        scores.append(score)

    sortScores = list(sorted(scores))
    cutoff = sortScores[-100]
    #cutoff = np.percentile(scores, perc)
    print "cutoff", cutoff

    filtSeqs = []
    for seq, score in zip(seqs, scores):
        if score > cutoff:
            filtSeqs.append(seq)
    print "sequences", len(filtSeqs)
            
    return filtSeqs

def makeLogo(seqs, title, fname):
    ofh = open("/tmp/seqs.fa", "w")
    for i, seq in enumerate(seqs):
        ofh.write(">%d\n%s\n" % (i, seq))
    ofh.close()
    print ofh.name

    cmd = "weblogo -t '%s' < /tmp/seqs.fa -F png -U probability > %s" % (title, fname)
    os.system(cmd)
    print "wrote %s" % fname

    os.remove("/tmp/seqs.fa")

def main():
    datasets = ["doench2014-Hs", "xu2015Train", "chari2015Train", "varshney2015", "gagnon2014"]
    fnames = []
    for dataset in datasets:
        seqs = parseTab("effData/%s.tab" % dataset, 75)
        fname = "out/logo-%s.png" % dataset
        makeLogo(seqs, dataset, fname)
        fnames.append(fname)

    cmd = "convert -append %s out/logos.png" % " ".join(fnames)
    os.system(cmd)


main()
