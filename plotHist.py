from annotateOffs import *
from collections import defaultdict
import operator

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.backends.backend_pdf as pltBack

def main():
    # parse offtargets.tsv
    offCounts = defaultdict(int)
    targetSeqs = dict()
    for row in iterTsvRows("offtargets.tsv"):
        offCounts[row.name]+=1

    xVals = offCounts.values()
    
    fname = "otCountHist.pdf"
    pdf = pltBack.PdfPages(fname)
    fig = plt.figure(figsize=(5,2),
                   dpi=300, facecolor='w')
    plt.hist(xVals, 80)
    fig.savefig(pdf, format = 'pdf')
    pdf.close()
    print fname

main()

