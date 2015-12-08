# a collection of off-target scoring functions from various websites/papers
# - MIT off-target score
# - CropIT
# - CFD
# - CCTop
# the central function is otScores()

# all scores have the property that the higher the score, the more likely is cutting
# one score (cctop) had to be inverted to follow this rule

import pickle, re
from os.path import basename, join, splitext, isfile, dirname

def otScores(guideSeq, otSeqs):
    """
    calculate all scores and return a dict scoreName -> list of scores
    scoreNames are: mit, cropit, cfd, cctop
    """
    scoreDesc = [
        ("cropit", calcCropitScore),
        ("cctop", calcCcTopScore),
        ("cfd", calcCfdScore),
        ("mit", calcMitScore)
    ]

    scoreDict = {}
    for scoreName, scoreFunc in scoreDesc:
        scores = []
        for otSeq in otSeqs:
            scores.append( scoreFunc(guideSeq, otSeq) )
        scoreDict[scoreName] = scores

    return scoreDict

def findRuns(lst):
    """ yield (start, end) tuples for all runs of ident. numbers in lst 
    >>> list(findRuns([1,1,1,0,0,1,0,1,1,1]))
    [(0, 3), (5, 6), (7, 10)]
    """
    start,end=False,False

    for i,x in enumerate(lst):
        if x and start is False:
            start=i
        if x==0 and start is not False and end is False:
            end=i-1
        if start is not False and end is not False:
            yield start,end+1       #and len is (end-start)
            start,end=False,False
    
    if start is not False:
        yield start,i+1       #and len is (end-start)

def calcCropitScore(guideSeq, otSeq):
    """
    see http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4605288/ PMID 26032770

    >>> int(calcCropitScore("GGGGGGGGGGGGGGGGGGGG","GGGGGGGGGGGGGGGGGGGG"))
    650

    # mismatch in 3' part
    >>> int(calcCropitScore("GGGGGGGGGGGGGGGGGGGA","GGGGGGGGGGGGGGGGGGGG"))
    575

    # mismatch in 5' part
    >>> int(calcCropitScore("AGGGGGGGGGGGGGGGGGGG","GGGGGGGGGGGGGGGGGGGG"))
    642

    # only mismatches -> least likely offtarget
    >>> int(calcCropitScore("AAAAAAAAAAAAAAAAAAAA","GGGGGGGGGGGGGGGGGGGG"))
    -27
    """
    if len(guideSeq)==23:
        guideSeq = guideSeq[:20]
        otSeq = otSeq[:20]

    assert(len(guideSeq)==len(otSeq)==20)

    penalties = [5,5,5,5,5,5,5,5,5,5,70,70,70,70,70,50,50,50,50,50]
    score = 0.0

    # do the score only for the non-mism positions
    misList = []
    score = 0.0
    for i in range(0, 20):
        if guideSeq[i]!=otSeq[i]:
            misList.append(1)
        else:
            misList.append(0)
            score += penalties[i]
    
    # get the runs of mismatches and update score for these positions
    consecPos = set()
    singlePos = set()
    for start, end in findRuns(misList):
        if end-start==1:
            score += -penalties[start]/2.0
        else:
            # mean if they happen to fall into different segments
            startScore = penalties[start]
            endScore = penalties[end-1]
            score += -((startScore+endScore)/2.0)

    return score

def calcCcTopScore(guideSeq, otSeq):
    """
    calculate the CC top score
    see http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0124633#sec002
    # no mismatch -> most likely off-target
    >>> int(calcCcTopScore("GGGGGGGGGGGGGGGGGGGG","GGGGGGGGGGGGGGGGGGGG"))
    224

    # mismatch in 5' part
    >>> int(calcCcTopScore("AGGGGGGGGGGGGGGGGGGG","GGGGGGGGGGGGGGGGGGGG"))
    222

    # mismatch in 3' part
    >>> int(calcCcTopScore("GGGGGGGGGGGGGGGGGGGA","GGGGGGGGGGGGGGGGGGGG"))
    185

    # only mismatches -> least likely offtarget
    >>> int(calcCcTopScore("AAAAAAAAAAAAAAAAAAAA","GGGGGGGGGGGGGGGGGGGG"))
    0
    """
    if len(guideSeq)==23:
        guideSeq = guideSeq[:20]
        otSeq = otSeq[:20]

    if not (len(guideSeq)==len(otSeq)==20):
        raise Exception("Not 20bp long: %s %dbp<-> %s %dbp" % (guideSeq, len(guideSeq), otSeq, len(otSeq)))
    score = 0.0
    for i in range(0, 20):
        if guideSeq[i]!=otSeq[i]:
            score += 1.2**(i+1)
    return 224.0-score

hitScoreM = [0,0,0.014,0,0,0.395,0.317,0,0.389,0.079,0.445,0.508,0.613,0.851,0.732,0.828,0.615,0.804,0.685,0.583]

def calcMitScore(string1,string2, startPos=0):
    """
    The MIT off-target score
    see 'Scores of single hits' on http://crispr.mit.edu/about
    startPos can be used to feed sequences longer than 20bp into this function

    the most likely off-targets have a score of 100
    >>> int(calcMitScore("GGGGGGGGGGGGGGGGGGGG","GGGGGGGGGGGGGGGGGGGG"))
    100

    # mismatches in the first three positions have no effect
    >>> int(calcMitScore("AGGGGGGGGGGGGGGGGGGG","GGGGGGGGGGGGGGGGGGGG"))
    100

    # less likely off-targets have lower scores
    >>> int(calcMitScore("GGGGGGGGGGGGGGGGGGGG","GGGGGGGGGGGGGGGGGGGA"))
    41
    """
    # The Patrick Hsu weighting scheme
    #print string1, string2
    if len(string1)==len(string2)==23:
        string1 = string1[:20]
        string2 = string2[:20]

    assert(len(string1)==len(string2)==20)

    dists = [] # distances between mismatches, for part 2
    mmCount = 0 # number of mismatches, for part 3
    lastMmPos = None # position of last mismatch, used to calculate distance

    score1 = 1.0
    for pos in range(0, len(string1)):
        if string1[pos]!=string2[pos]:
            mmCount+=1
            if lastMmPos!=None:
                dists.append(pos-lastMmPos)
            score1 *= 1-hitScoreM[pos]
            lastMmPos = pos

    # 2nd part of the score - distribution of mismatches
    if mmCount<2: # special case, not shown in the paper
        score2 = 1.0
    else:
        avgDist = sum(dists)/len(dists)
        score2 = 1.0 / (((19-avgDist)/19.0) * 4 + 1)

    # 3rd part of the score - mismatch penalty
    if mmCount==0: # special case, not shown in the paper
        score3 = 1.0
    else:
        score3 = 1.0 / (mmCount**2)

    score = score1 * score2 * score3 * 100
    return score

compTable = { "a":"t", "A":"T", "t" :"a", "T":"A", "c":"g", "C":"G", "g":"c", "G":"C", "N":"N", "n":"n", 
        "Y":"R", "R" : "Y", "M" : "K", "K" : "M", "W":"W", "S":"S",
        "H":"D", "B":"V", "V":"B", "D":"H", "y":"r", "r":"y","m":"k",
        "k":"m","w":"w","s":"s","h":"d","b":"v","d":"h","v":"b","y":"r","r":"y" }

def complRna(seq):
    " complement the sequence and translate to RNA "
    newseq = []
    for nucl in seq.upper():
        newseq.append( compTable[nucl].replace("T", "U") )
    return "".join(newseq)

def outMats():
    """
    write normalized matrices to out 
    """
    prettyMat("all", "out/hsuAll.tsv")
    prettyMat("row", "out/hsuRow.tsv")
    prettyMat("col", "out/hsuCol.tsv")
    prettyMat("none", "out/hsuNone.tsv")

def prettyMat(strat, outFname):
    """
    wrote normalized Hsu matrix to outFname
    """
    ofh = open(outFname, "w")
    row = ["pos"]
    row.extend(range(2,21))
    row = [str(x) for x in row]
    ofh.write( "\t".join(row)+"\n")

    normMat, normAvg  = parseHsuMat("./hsu2013/fig2cData.txt", strat=strat)
    for nucl, freqs in normMat.iteritems():
        row = [":".join(nucl)]
        row.extend(freqs)
        row = [str(x) for x in row]
        ofh.write( "\t".join(row)+"\n")

    row = ["avg"]
    row.extend(normAvg)
    row = [str(x) for x in row]
    ofh.write( "\t".join(row)+ "\n")

rawHsuMat= {('A', 'A'): [0.271959618, 0.697036065, 0.171770967,
0.021754414, 0.005241839, 0.025188071, 0.046859394, 0.267506677,
0.259177567, 0.07797037, 1.286750463, 0.535596183, 1.032247542,
0.643922558, 0.515958107, 1.130955312, 2.870906184, 1.095411389,
1.191602132], ('G', 'G'): [0.496472865, 0.525028618, 0.215301755,
0.226670927, 0.004546068, 0.04966372, 0.089927088, 0.348976795,
1.579071836, 0.707528464, 1.155940112, 0.810035714, 1.016250559,
0.230357513, 0.95239301, 1.205954808, 1.296737018, 1.415685748,
1.323358326], ('U', 'T'): [0.337223457, 0.168000344, 0.115211807,
0.261310401, 0.279295225, 0.266509558, 0.179895816, 0.619215107,
0.356050082, 1.14001417, 0.62100718, 0.299685354, 0.836770787,
0.900417256, 0.738011252, 1.181552961, 1.157934409, 0.821526176,
1.734698905], ('U', 'C'): [0.254333541, 0.085714515, 0.054631234,
0.36587133, 0.016462872, 0.068822052, 0.046215969, 0.09350449,
0.263954537, 0.340470548, 0.214279635, 0.436699042, 2.340127222,
0.54281779, 0.010174614, 1.341048155, 1.646620387, 0.995399409,
1.428913939], ('C', 'T'): [0.440303113, 0.135104168, 0.043264802,
0.242026012, 0.327875553, 0.072221711, 0.137942313, 0.5575506,
0.432695002, 0.904623856, 0.448727403, 0.305879186, 0.863503533,
1.014157582, 0.647532031, 1.086889349, 1.073187004, 0.968782834,
1.671881765], ('G', 'A'): [0.458774127, 0.562057357, 0.317460049,
0.368206347, 0.047756641, 0.022744096, 0.01322757, 0.249708488,
0.195339521, 0.123994188, 1.365452127, 0.688216641, 0.787503713,
0.813153469, 0.60357013, 1.046113038, 2.786745287, 0.454444131,
1.16622155], ('G', 'T'): [0.689780687, 0.282072211, 0.537606185,
0.839891532, 0.716358798, 0.898073668, 0.52549235, 1.043403879,
0.694185851, 1.521066692, 0.711884156, 0.26584025, 1.01625944,
0.963232124, 0.807022991, 0.923530678, 1.597721072, 0.903020898,
1.653312217], ('C', 'A'): [0.166382281, 0.531452139, 0.178570178,
0.414700936, 0.50126924, 0.036677831, 0.116486245, 0.726931853,
0.617817481, 0.35333324, 1.4293465, 0.110193397, 1.065079759,
0.75614004, 0.240529657, 1.622197865, 2.250415945, 0.964947064,
1.018630984], ('A', 'G'): [0.027493368, 0.217519417, 0.100600517,
0.155294081, 0.025303057, 0.059154973, 0.036503049, 0.03030029,
0.820965614, 0.304721275, 1.037514324, 0.602176238, 0.927208246,
0.491187717, 0.976263629, 0.744956324, 1.295205389, 1.082002774,
1.478634617], ('U', 'G'): [0.060056479, 0.134281319, 0.311857733,
0.226053995, 0.054924206, 0.071959174, 0.128509362, 0.416561912,
1.73814687, 0.718418143, 1.53185795, 0.988456707, 1.24449999,
0.779052105, 1.064452173, 1.219746615, 1.268864823, 1.27712259,
1.536590735], ('C', 'C'): [0.059633277, 0.025314601, 0.055185616,
0.211364672, 0.009880719, 0.045811702, 0.010531506, 0.041691214,
0.10261082, 0.136908112, 0.155528383, 0.458455298, 2.080996466,
0.639614156, 0.054046079, 1.407539975, 1.050709161, 1.047816491,
1.428913939], ('A', 'C'): [0.270019825, 0.136179168, 0.181070971,
0.349616508, 0.018107063, 1.022086513, 0.705976348, 0.482360483,
0.48415093, 0.729154118, 0.927047872, 0.773675684, 1.8686911,
0.552164002, 0.024710323, 0.992387533, 1.658152541, 1.159682815,
1.514208159]}

def calcRawHsu(guideSeq, otSeq):
    """ raw sum of hsu freqs 
    >>> calcRawHsu("GAGTCCGAGCAGAAGAAGAA","GAGTCAGAACAGAAGAACAA")
    0.0
    >>> calcRawHsu("GGGGGGGGGGGGGGGGGGGG","GGGGGGGGGGGGGGGGGGGG")
    0.0
    >>> calcRawHsu("GGGGGGGGGGGGGGGGGGGG","GGGGGGGGGGGGGGGGGGGT")
    0.0
    >>> calcRawHsu("GGGGGGGGGGGGGGGGGGGG","GGGGGGGGGGGGGGGGGGGC")
    0.0
    >>> calcRawHsu("GGGGGGGGGGGGGGGGGGGG","GGGGGGGGGGGGGGGGGCGG")
    0.0

    """
    guideSeq = guideSeq[1:20]
    otSeq = otSeq[1:20]
    #print guideSeq, len(guideSeq), otSeq, len(otSeq)
    assert(len(guideSeq)==19)
    assert(len(otSeq)==19)# Hsu ignores pos 0
    rnaSeq = complRna(guideSeq)

    mismatchPosList = []
    score = 1.0
    for i in range(0, 19):
        rnaNucl, dnaNucl = rnaSeq[i], otSeq[i]
        # "In case of a match, both were set equal to 1."
        if (rnaNucl, dnaNucl) in [('C', 'G'), ('U', 'A'), ('A', 'T'), ('G', 'C')]:
            f = 1.0
        else:
            f = rawHsuMat[(rnaNucl, dnaNucl)][i]
        score *= f
    return score

hsuMat = None # dict with (fromNucl, toNucl) -> list of 19 scores
avgFreqs = None # list of 19 scores
hsuStrat = None # loaded matrix has strat

# see hsuMat.py
nuclFreqs= {('A', 'A'): 0.4819440141370613, ('G', 'G'): 0.6571297543038187, ('U', 'T'): 0.4663759334374003, ('U', 'C'): 0.18755561795635842, ('C', 'T'): 0.3917125484856841, ('G', 'A'): 0.472948896301865, ('G', 'T'): 1.0, ('A', 'G'): 0.2796160896968995, ('U', 'G'): 0.787929779020387, ('C', 'C'): 0.0, ('A', 'C'): 0.6804018833297995, ('C', 'A'): 0.5931243444910334}
posFreqs = [0.294369386, 0.29164666, 0.190210984, 0.306896763, 0.167251773, 0.219909422, 0.169797251, 0.406475982, 0.628680509, 0.588183598, 0.907111342, 0.522909141, 1.256594863, 0.693851359, 0.552888666, 1.158572718, 1.662766602, 1.01548686, 1.428913939]

def calcHsuSuppScore2(guideSeq, otSeq):
    """
    The Hsu score in a version that only uses the aggregrate frequencies

    >>> calcHsuSuppScore2("GAGTCCGAGCAGAAGAAGAA","GAGTCCGAGCAGAAGAAGAG")
    1.1258838441954209
    >>> calcHsuSuppScore2("GTGTCCGAGCAGAAGAAGAA","GAGTCCGAGCAGAAGAAGAA")
    0.14186956352790206
    >>> calcHsuSuppScore2("GAGTCCGAGCAGAAGAAGAA","GAGTCAGAACAGAAGAACAA")
    0.0
    >>> calcHsuSuppScore2("GGGGGGGGGGGGGGGGGGGG","GGGGGGGGGGGGGGGGGGGG")
    1.0
    >>> calcHsuSuppScore2("GGGGGGGGGGGGGGGGGGGG","GGGGGGGGGGGGGGGGGGGT")
    0.5597235206124074
    >>> calcHsuSuppScore2("GGGGGGGGGGGGGGGGGGGG","GGGGGGGGGGGGGGGGGGGC")
    0.0
    >>> calcHsuSuppScore2("GGGGGGGGGGGGGGGGGGGG","GGGGGGGGGGGGGGGGGCGG")
    0.0
    >>> calcHsuSuppScore2("GAGGGGGGGGGGGGGGGGGG","GGGGGGGGGGGGGGGGGGGG")
    0.231942405261347
    """
    guideSeq = guideSeq[1:20]
    otSeq = otSeq[1:20]
    #print guideSeq, len(guideSeq), otSeq, len(otSeq)
    assert(len(guideSeq)==19)
    assert(len(otSeq)==19)# Hsu ignores pos 0

    rnaSeq = complRna(guideSeq)
    # "Predicted cutting frequencies for genome-wide targets were calculated by
    # multiplying, in series: fest = f(1) * g(N1,N1') * f(2) * g(N2,N2') * ... * h
    # with values f(i) and g(Ni, Ni') at position i corresponding,
    # respectively, to the aggregate position- and base-mismatch cutting
    # frequencies for positions and pairings indicated in Fig. 2c"
    mismatchPosList = []
    score = 1.0
    for i in range(0, 19):
        rnaNucl, dnaNucl = rnaSeq[i], otSeq[i]
        # "In case of a match, both were set equal to 1."
        if (rnaNucl, dnaNucl) in [('C', 'G'), ('U', 'A'), ('A', 'T'), ('G', 'C')]:
            f = 1.0
            g = 1.0
        else:
            f = posFreqs[i]
            g = nuclFreqs[(rnaNucl, dnaNucl)]
            mismatchPosList.append(i)
        score *= f * g

    # "the value h meanwhile re-weighted the estimated
    # frequency by the minimum pairwise distance between consecutive mismatches
    # in the target sequence. this distance value, in base-pairs, was divided
    # by 18 to give a maximum value of 1 (in cases where fewer than 2
    # mismatches existed, or where mismatches occurred on opposite ends of
    # the 19 bp target-window"
    if len(mismatchPosList)<2:
        h = 1.0
    else:
        dists = []
        for left, right in zip(mismatchPosList[:-1], mismatchPosList[1:]):
            dists.append(right-left)
        minDist = min(dists)
        h = minDist / 18.0
    score *= h

    return score

def loadHsuMat(strat):
    if strat in ["avgs", "raw"]:
        return
    global hsuMat
    global avgFreqs
    #global hsuStrat
    matFname = "./hsu2013/fig2cData.txt"
    hsuMat, avgFreqs = parseHsuMat(matFname, strat)
    #hsuStrat = strat

def parseHsuMat(fname, strat="col"):
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

def calcHsuSuppScore(guideSeq, otSeq, baseDir="./", strat="all", dfh=None):
    """ calculate the score described on page 17 of the Hsu et al 2013 supplement PDF
    Will ignore position 0 of both the ot and the guide as the Hsu score is only
    defined for positions 1-20

    # mismatch in 5' part -> 

    >> calcHsuSuppScore("GAGTCCGAGCAGAAGAAGAA","GAGTCCGAGCAGAAGAAGAG")
    0.4509132855355929
    >> calcHsuSuppScore("GTGTCCGAGCAGAAGAAGAA","GAGTCCGAGCAGAAGAAGAA")
    0.007929899123452079
    >> calcHsuSuppScore("GAGTCCGAGCAGAAGAAGAA","GAGTCAGAACAGAAGAACAA")
    3.69683024017458e-08
    >> calcHsuSuppScore("GGGGGGGGGGGGGGGGGGGG","GGGGGGGGGGGGGGGGGGGG")
    1.0
    >> calcHsuSuppScore("GGGGGGGGGGGGGGGGGGGG","GGGGGGGGGGGGGGGGGGGT")
    0.4907323091938597
    >> calcHsuSuppScore("GGGGGGGGGGGGGGGGGGGG","GGGGGGGGGGGGGGGGGGGC")
    0.00880386936389443
    >> calcHsuSuppScore("GAGGGGGGGGGGGGGGGGGG","GGGGGGGGGGGGGGGGGGGG")
    0.0016461094044462768
    """
    if strat=="avgs":
        return calcHsuSuppScore2(guideSeq, otSeq)
    if strat=="raw":
        return calcRawHsu(guideSeq, otSeq)

    guideSeq = guideSeq[1:20]
    otSeq = otSeq[1:20]
    #print guideSeq, len(guideSeq), otSeq, len(otSeq)
    assert(len(guideSeq)==19)
    assert(len(otSeq)==19)# Hsu ignores pos 0

    global hsuMat
    global avgFreqs
    global hsuStrat
    if hsuMat is None or hsuStrat!=strat:
        matFname = baseDir+"hsu2013/fig2cData.txt"
        hsuMat, avgFreqs = parseHsuMat(matFname, strat)
        hsuStrat = strat

    rnaSeq = complRna(guideSeq)
    if dfh:
        dfh.write("guideDna=%s, guideRna=%s, otSeq=%s\n" % (guideSeq, rnaSeq, otSeq))
    # "Predicted cutting frequencies for genome-wide targets were calculated by
    # multiplying, in series: fest = f(1) * g(N1,N1') * f(2) * g(N2,N2') * ... * h
    # with values f(i) and g(Ni, Ni')
    # at position i corresponding, respectively, to the aggregate
    # position- and base-mismatch cutting frequencies for positions and pairings indicated in Fig. 2c"
    mismatchPosList = []
    score = 1.0
    for i in range(0, 19):
        rnaNucl, dnaNucl = rnaSeq[i], otSeq[i]
        # "In case of a match, both were set equal to 1."
        if (rnaNucl, dnaNucl) in [('C', 'G'), ('U', 'A'), ('A', 'T'), ('G', 'C')]:
            f = 1.0
            g = 1.0
        else:
            f = avgFreqs[i]
            g = hsuMat[(rnaNucl, dnaNucl)][i]
            mismatchPosList.append(i)
        if strat.endswith("_Sum"):
            score += f*g
        elif strat.endswith("_allSum"):
            score += f+g
        else:
            score *= f * g
        if dfh:
            dfh.write("pos: %d, RNA-nucl: %s, DNA-nucl: %s -> f=%f, g=%f, score=%f\n" % ( i, rnaNucl, dnaNucl, f, g, score))
        #score *= g

    # "The value h meanwhile re-weighted the estimated
    # frequency by the minimum pairwise distance between consecutive mismatches in the target
    # sequence. This distance value, in base-pairs, was divided by 18 to give a maximum value of 1 (in
    # cases where fewer than 2 mismatches existed, or where mismatches occurred on opposite ends of
    # the 19 bp target-window"
    if len(mismatchPosList)<2:
        h = 1.0
    else:
        dists = []
        for left, right in zip(mismatchPosList[:-1], mismatchPosList[1:]):
            dists.append(right-left)
        minDist = min(dists)
        h = minDist / 18.0
    score *= h

    if dfh:
        dfh.write("h=%f\n" % h)
        dfh.write("score=%f\n\n" % score)
    return score

# === SOURCE CODE cfd-score-calculator.py provided by John Doench =====

def get_mm_pam_scores():
    """
    """
    dataDir = join(dirname(__file__), 'CFD_Scoring')
    mm_scores = pickle.load(open(join(dataDir, 'mismatch_score.pkl'),'rb'))
    pam_scores = pickle.load(open(join(dataDir, 'pam_scores.pkl'),'rb'))
    return (mm_scores,pam_scores)
    #except:
        #raise Exception("Could not find file with mismatch scores or PAM scores")

#Reverse complements a given string
def revcom(s):
    basecomp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','U':'A'}
    letters = list(s[::-1])
    letters = [basecomp[base] for base in letters]
    return ''.join(letters)

#Calculates CFD score
def calc_cfd(wt,sg,pam):
    """
    >>> calc_cfd("GGGGGGGGGGGGGGGGGGGG", "GGGGGGGGGGGGGGGGGGGG", "GG")
    1.0
    """
    mm_scores,pam_scores = get_mm_pam_scores()
    score = 1
    sg = sg.replace('T','U')
    wt = wt.replace('T','U')
    s_list = list(sg)
    wt_list = list(wt)
    #print mm_scores
    for i,sl in enumerate(s_list):
        #print i, sl, wt_list[i], mm_scores
        if wt_list[i] == sl:
            score*=1
        else:
            key = 'r'+wt_list[i]+':d'+revcom(sl)+','+str(i+1)
            score *= mm_scores[key]
    score*=pam_scores[pam]
    return (score)

mm_scores, pam_scores = None, None

def calcCfdScore(guideSeq, otSeq):
    """ based on source code provided by John Doench
    >>> calcCfdScore("GGGGGGGGGGGGGGGGGGGGGGG", "GGGGGGGGGGGGGGGGGAAAGGG")
    0.4635989007074176
    >>> calcCfdScore("GGGGGGGGGGGGGGGGGGGGGGG", "GGGGGGGGGGGGGGGGGGGGGGG")
    1.0
    >>> calcCfdScore("GGGGGGGGGGGGGGGGGGGGGGG", "aaaaGaGaGGGGGGGGGGGGGGG")
    0.5140384614450001
    """
    global mm_scores, pam_scores
    if mm_scores is None:
        mm_scores,pam_scores = get_mm_pam_scores()
    wt = guideSeq.upper()
    off = otSeq.upper()
    m_wt = re.search('[^ATCG]',wt)
    m_off = re.search('[^ATCG]',off)
    if (m_wt is None) and (m_off is None):
        pam = off[-2:]
        sg = off[:-3]
        cfd_score = calc_cfd(wt,sg,pam)
        return cfd_score
        #print "CFD score: "+str(cfd_score)

# ==== END CFD score source provided by John Doench

#outMats()
#dfh = open("out/hsuDebug.txt", "w")
#calcHsuSuppScore("GGGGGGGGGGGGGGGGGGGG","GGGGGGGGGGGGGGGGGGGG", dfh=dfh)
#calcHsuSuppScore("GAGTCCGAGCAGAAGAAGAA","GAGTCCGAGCAGAAGAAGAG", dfh=dfh)
#calcHsuSuppScore("GTGTCCGAGCAGAAGAAGAA","GAGTCCGAGCAGAAGAAGAA", dfh=dfh)
#calcHsuSuppScore("GTGTCCGAGCAGAAGAAGAA","GGGGGGGGGGGGGGGGGGGG", dfh=dfh)
#loadHsuMat("none")
#print avgFreqs
if __name__=="__main__":
    import doctest
    doctest.testmod()
