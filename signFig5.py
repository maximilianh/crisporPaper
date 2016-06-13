# calculate if the differences from fig5 are significant

import math, os
import numpy as np
import annotateOffs
from scipy import stats
#import rpy
#import rpy2.robjects as ro

fnames = ["effData/concordet2.scores.tab", "effData/alenaAll.scores.tab", \
    "effData/teboulVivo_mm9.scores.tab", "effData/schoenig.scores.tab", "effData/eschstruth.scores.tab"]

#corrOverlap = ro.r['cocor.dep.groups.overlap']

ofh = open("out/signFig5.tab", "w")

for fname in fnames:
        x1 = []
        x2 = []
        y = []
        #for scoreType in ["fusi", "crisprScan"]:
        for row in annotateOffs.iterTsvRows(fname):
            x1.append(float(row.fusi))
            x2.append(float(row.crisprScan))
            y.append(float(row.modFreq))

        r1 = stats.pearsonr(x1, y)[0]
        r2 = stats.pearsonr(x2, y)[0]
        r12 = stats.pearsonr(x1, x2)[0]
        #res = corrOverlap(r.jk=r1, r.jh=r2, r.kh=r12, n=len(x1), alternative="less", alpha=0.05, conf.level=0.95, null.value=0)
        if r1 < r2:
            way = "less"
        else:
            way = "more"
        row = [fname, r1, r2, r12, len(x1), way]
        row = [str(x) for x in row]
        ofh.write("\t".join(row)+"\n")

        # confidence intervals
        #num = float(len(x))
        #stderr = 1.0 / math.sqrt(num - 3)
        #delta = 1.96 * stderr
        #lower = math.tanh(math.atanh(r) - delta)
        #upper = math.tanh(math.atanh(r) + delta)
        #print fname, len(x1), r1, r2, r3

    #z = np.arctanh(corr[0])
    #std = np.std(x)

os.system("Rscript signFig5calc.R")


