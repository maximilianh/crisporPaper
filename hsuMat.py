# output the matrices from hsu 2013 in a format suitable for pasting into 
# python code

# from suppl file S2, Excel sheet 1

def parseHsuMat(fname, strat="all"):
    """ return the hsu 2013 matrix as a dict rnaNucl -> dnaNucl -> list of scores and a list of 19 averages, one per position
    >> parseHsuMat("./hsu2013/fig2cData.txt")
    """
    hsuMat = {}
    minMat = 99999.0
    maxMat = 0.0
    for line in open(fname):
        if line.startswith("nucl"):
            continue
        fs = line.rstrip("\n").split()
        # the values are in the order 19-1 3'-5' in the file, but our sequences are always 1-19, 5'-3'
        freqs = list(reversed([float(x) for x in fs[1:]]))
        if line.startswith("avg"):
            avgs = freqs
            continue
        nuclComb = fs[0]
        rnaNucl, dnaNucl = nuclComb.split(":")
        hsuMat[ (rnaNucl, dnaNucl)] = freqs
        minMat = min(min(freqs), minMat)
        maxMat = max(max(freqs), maxMat)

    minCols = []
    maxCols = []
    for i in range(0, 19):
        colVals = []
        for freqs in hsuMat.values():
            colVals.append(freqs[i])
        minCols.append(min(colVals))
        maxCols.append(max(colVals))

    pCount = 0.0001 # must use pseudo counts
    # normalize
    # "Each frequency was normalized to range from 0 to 1, such that f = (f-fmin) / (fmax-fmin)"
    normMat = {}
    for key, freqs in hsuMat.iteritems():
        if strat=="all":
            normVals = [( f - minMat) / (maxMat - minMat) for f in freqs]
        elif strat=="row":
            minFreq, maxFreq = min(freqs), max(freqs)
            normVals = [( f - minFreq) / (maxFreq - minFreq) for f in freqs]
        elif strat=="col":
            normVals = [( f - minCols[i]) / (maxCols[i] - minCols[i]) for i, f in enumerate(freqs)]
        elif strat.startswith("none") or strat=="onlyAvgs":
            normVals = freqs
        elif strat.startswith("limit"):
            normVals = [min(f, 1.0) for f in freqs]
        else:
            assert(False)

        if not strat.startswith("none"):
            normVals = [pCount+n for n in normVals]
        normMat[key] = normVals
    hsuMat = normMat

    if strat in ["all", "row", "onlyAvgs", "limit"]:
        normAvgs = [(a-min(avgs)) / (max(avgs)-min(avgs)) for a in avgs]
    else:
        normAvgs = avgs

    if not strat.startswith("none"):
        normAvgs = [n+pCount for n in normAvgs]
    assert(min(normAvgs)!=0.0)
    avgs = normAvgs

    #print "loaded hsu mat", strat
    #print hsuMat
    #print avgs
    return hsuMat, avgs

def normalizeArr(arr):
    """
    normalize a list of numbers as in Hsu 2013
    """

# last row
posFreqs = [0.294369386, 0.29164666, 0.190210984, 0.306896763, 0.167251773, 0.219909422, 0.169797251, 0.406475982, 0.628680509, 0.588183598, 0.907111342, 0.522909141, 1.256594863, 0.693851359, 0.552888666, 1.158572718, 1.662766602, 1.01548686, 1.428913939]
minX = min(posFreqs)
maxX = max(posFreqs)
arrDiff = maxX-minX
posFreqs = [ x - minX / arrDiff for x in posFreqs]

# last column
nuclFreqs = {
    ("G","T") : 0.873145036,
    ("A","C") : 0.728917998,
    ("G","G") : 0.718415839,
    ("U","G") : 0.777442783,
    ("U","T") : 0.632333171,
    ("C","A") : 0.689531718,
    ("C","T") : 0.598639359,
    ("G","A") : 0.635299393,
    ("A","A") : 0.639358676,
    ("A","G") : 0.54805289,
    ("U","C") : 0.506508186,
    ("C","C") : 0.421868791,
}
minFreq = min(nuclFreqs.values())
maxFreq = max(nuclFreqs.values())
maxDiff = maxFreq-minFreq

newFreqs = {}
for nucl, freq in nuclFreqs.iteritems():
    newFreqs[nucl] = (freq-minFreq)/maxDiff

hsuMat, avgs = parseHsuMat("./hsu2013/fig2cData.txt", strat="none")

#print "posFreqs=",posFreqs
#print "nuclFreqs=",newFreqs
print "hsuMat=", hsuMat
