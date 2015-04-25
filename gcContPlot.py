from annotateOffs import *
from collections import defaultdict
import operator

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.backends.backend_pdf as pltBack

def plot(fname):
    outfname = "gcCont" + '.pdf'
    pdf = pltBack.PdfPages(outfname)
    fig = plt.figure(figsize=(5,5),
                   dpi=300, facecolor='w')
    data = np.genfromtxt(fname, names=True, dtype=np.dtype([('name', 'S20'), ('gcContent', 'f8'), ('offtargetCount', 'f8')]))
    fig = plt.figure()

    studyX = defaultdict(list)
    studyY = defaultdict(list)
    for row in data:
        study = row["guideName"].split("_")[0]
        gcCont = row["gcContent"]
        otCount = row["offtargetCount"]
        studyX[study].append( gcCont )
        studyY[study].append( otCount )
    
    colors = ["green", "blue", "black", "yellow", "red", "grey"]
    markers = ["o", "s", "+", ">", "<", "^"]
    studyNames = []
    figs = []
    i = 0
    for study, xVals in studyX.iteritems():
        yVals = studyY[study]
        studyFig = plt.scatter(xVals, yVals, \
            alpha=.5, \
            marker=markers[i], \
            s=60, \
            color=colors[i])
        figs.append(studyFig)
        studyNames.append(study)
        i+=1
    
    plt.legend(figs,
           studyNames,
           scatterpoints=1,
           loc='upper left',
           ncol=3,
           fontsize=10)

    plt.xlabel("GC content")
    plt.ylabel("Number of off-targets")
    fig.savefig(pdf, format = 'pdf')
    pdf.close()
    print "Wrote %s" % outfname



def main():
    # parse offtargets.tsv
    offsByName = defaultdict(list)
    targetSeqs = dict()
    for row in iterTsvRows("offtargets.tsv"):
        offsByName[row.name].append(row)
        if row.type=="on-target":
            targetSeqs[row.name] = row.seq

    rows = []
    for name, offs in offsByName.iteritems():
        row = [name, str(gcCont(targetSeqs[name])), str(len(offs))]
        rows.append(row)
    rows.sort(key=operator.itemgetter(1))

    # write to tsv file
    fname = "gcCont.tsv"
    ofh = open(fname, "w")
    headers = ["guideName", "gcContent", "offtargetCount"]
    ofh.write("\t".join(headers))
    ofh.write("\n")

    for row in rows:
        ofh.write("\t".join(row))
        ofh.write("\n")
    ofh.close()

    plot(fname)

main()
