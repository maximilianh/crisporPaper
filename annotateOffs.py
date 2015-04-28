# annotate guideseq offtargets with all possible scores we have
import glob, copy, sys, math, operator, random, re, collections
from collections import defaultdict, Counter
from os.path import basename, join

def iterTsvRows(inFile, encoding=None, fieldSep="\t", isGzip=False, skipLines=None, \
        makeHeadersUnique=False, commentPrefix=None, headers=None):
    """ 
        parses tab-sep file with headers as field names 
        yields collection.namedtuples
        strips "#"-prefix from header line
    """

    if isinstance(inFile, str):
        if inFile.endswith(".gz") or isGzip:
            fh = gzip.open(inFile, 'rb')
        else:
            fh = open(inFile)
    else:
        fh = inFile

    if headers==None:
        line1 = fh.readline()
        line1 = line1.strip("\n").strip("#")
        headers = line1.split(fieldSep)
        headers = [re.sub("[^a-zA-Z0-9_]","_", h) for h in headers]

    if makeHeadersUnique:
        newHeaders = []
        headerNum = defaultdict(int)
        for h in headers:
            headerNum[h]+=1
            if headerNum[h]!=1:
                h = h+"_"+str(headerNum[h])
            newHeaders.append(h)
        headers = newHeaders

    if skipLines:
        for i in range(0, skipLines):
            fh.readline()

    Record = collections.namedtuple('tsvRec', headers)
    for line in fh:
        if commentPrefix!=None and line.startswith(commentPrefix):
            continue
        line = line.rstrip("\n")
        fields = line.split(fieldSep)
        if encoding!=None:
            fields = [f.decode(encoding) for f in fields]
        #fields = [x.decode(encoding) for x in fields]
        try:
            rec = Record(*fields)
        except Exception, msg:
            logging.error("Exception occured while parsing line, %s" % msg)
            logging.error("Filename %s" % fh.name)
            logging.error("Line was: %s" % line)
            logging.error("Does number of fields match headers?")
            logging.error("Headers are: %s" % headers)
            raise Exception("wrong field count in line %s" % line)
        # convert fields to correct data type
        yield rec

def parseGuides():
    " return guides as dict name -> seq "
    guides = {}
    guidesExt = {}
    for l in open("guides.txt"):
        name, seq, extSeq, pos = l.strip().split()
        guides[name] = seq
        guidesExt[name] = extSeq
    return guides, guidesExt
        
hitScoreM = [0,0,0.014,0,0,0.395,0.317,0,0.389,0.079,0.445,0.508,0.613,0.851,0.732,0.828,0.615,0.804,0.685,0.583]

def calcHitScore(string1,string2, startPos=0):
    """ see 'Scores of single hits' on http://crispr.mit.edu/about 
    startPos can be used to feed sequences longer than 20bp into this function
    """
    # The Patrick Hsu weighting scheme
    #print string1, string2
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
    # 2nd part of the score
    if mmCount<2: # special case, not shown in the paper
        score2 = 1.0
    else:
        avgDist = sum(dists)/len(dists)
        score2 = 1.0 / (((19-avgDist)/19.0) * 4 + 1)
    # 3rd part of the score
    if mmCount==0: # special case, not shown in the paper
        score3 = 1.0
    else:
        score3 = 1.0 / (mmCount**2)

    score = score1 * score2 * score3 * 100
    return score


def parseBeds(dirName):
    " parse all beds in dir and return as baseFname -> set of (seq,score, count,guideName)"
    ret = {}
    for fname in glob.glob(dirName+"/*.bed"):
        counts = defaultdict(int)
        scores = defaultdict(int)
        for line in open(fname):
            fs = line.rstrip("\n").split("\t")
            seq = fs[3]
            score = float(fs[4])
            counts[seq]+=1
            scores[seq]+=score

        guideName = basename(fname).split('.')[0]
        seqs = set()
        for seq, count in counts.iteritems():
            score = scores[seq]
            seqs.add( (seq, score, count, guideName) )

        ret[guideName]= seqs
    return ret

def parseFasta(fileObj):
    " parse a fasta file, where each seq is on a single line, return dict id -> seq "
    seqs = {}
    parts = []
    seqId = None
    for line in fileObj:
        line = line.rstrip("\n")
        if line.startswith(">"):
            if seqId!=None:
                seqs[seqId]  = "".join(parts)
            seqId = line.lstrip(">")
            parts = []
        else:
            parts.append(line)
    if len(parts)!=0:
        seqs[seqId]  = "".join(parts)
    return seqs

def parseFlanks(faDir):
    seqToFlank = {}
    for fname in glob.glob(join(faDir, "*.fa")):
        for shortSeq, seq in parseFasta(open(fname)).iteritems():
            assert(shortSeq not in seqToFlank)
            assert(len(seq)==100)
            seqToFlank[shortSeq] = seq
    return seqToFlank

def gcCont(seq):
    seq = seq.upper()
    return int(100*float(seq.count("C")+seq.count("G")) / len(seq))

# DOENCH SCORING 
params = [
# pasted/typed table from PDF and converted to zero-based positions
(1,'G',-0.2753771),(2,'A',-0.3238875),(2,'C',0.17212887),(3,'C',-0.1006662),
(4,'C',-0.2018029),(4,'G',0.24595663),(5,'A',0.03644004),(5,'C',0.09837684),
(6,'C',-0.7411813),(6,'G',-0.3932644),(11,'A',-0.466099),(14,'A',0.08537695),
(14,'C',-0.013814),(15,'A',0.27262051),(15,'C',-0.1190226),(15,'T',-0.2859442),
(16,'A',0.09745459),(16,'G',-0.1755462),(17,'C',-0.3457955),(17,'G',-0.6780964),
(18,'A',0.22508903),(18,'C',-0.5077941),(19,'G',-0.4173736),(19,'T',-0.054307),
(20,'G',0.37989937),(20,'T',-0.0907126),(21,'C',0.05782332),(21,'T',-0.5305673),
(22,'T',-0.8770074),(23,'C',-0.8762358),(23,'G',0.27891626),(23,'T',-0.4031022),
(24,'A',-0.0773007),(24,'C',0.28793562),(24,'T',-0.2216372),(27,'G',-0.6890167),
(27,'T',0.11787758),(28,'C',-0.1604453),(29,'G',0.38634258),(1,'GT',-0.6257787),
(4,'GC',0.30004332),(5,'AA',-0.8348362),(5,'TA',0.76062777),(6,'GG',-0.4908167),
(11,'GG',-1.5169074),(11,'TA',0.7092612),(11,'TC',0.49629861),(11,'TT',-0.5868739),
(12,'GG',-0.3345637),(13,'GA',0.76384993),(13,'GC',-0.5370252),(16,'TG',-0.7981461),
(18,'GG',-0.6668087),(18,'TC',0.35318325),(19,'CC',0.74807209),(19,'TG',-0.3672668),
(20,'AC',0.56820913),(20,'CG',0.32907207),(20,'GA',-0.8364568),(20,'GG',-0.7822076),
(21,'TC',-1.029693),(22,'CG',0.85619782),(22,'CT',-0.4632077),(23,'AA',-0.5794924),
(23,'AG',0.64907554),(24,'AG',-0.0773007),(24,'CG',0.28793562),(24,'TG',-0.2216372),
(26,'GT',0.11787758),(28,'GG',-0.69774)]

intercept =  0.59763615
gcHigh    = -0.1665878
gcLow     = -0.2026259

def calcDoenchScore(seq):
    assert(len(seq)==30)
    score = intercept

    guideSeq = seq[4:24]
    gcCount = guideSeq.count("G") + guideSeq.count("C")
    if gcCount <= 10:
        gcWeight = gcLow
    if gcCount > 10:
        gcWeight = gcHigh
    score += abs(10-gcCount)*gcWeight

    for pos, modelSeq, weight in params:
        subSeq = seq[pos:pos+len(modelSeq)]
        if subSeq==modelSeq:
            score += weight
    return 1.0/(1.0+math.exp(-score))

# Microhomology score from Bae et al, Nat Biotech 2014 

def calcMicroHomolScore(seq, left):
    """ calculate the micro homology and out-of-frame score for a breakpoint in a 60-80mer
    See http://www.nature.com/nmeth/journal/v11/n7/full/nmeth.3015.html
    Source code adapted from Supp File 1

    From the manuscript:
    "On the basis of these observations, we developed a simple formula and a
    computer program (Supplementary Fig. 3) to predict the deletion patterns
    at a given nuclease target site that are associated with microhomology of
    at least two bases (Fig. 1b and Supplementary Note). We assigned a pattern
    score to each deletion pattern and a microhomology score (equaling the sum
    of pattern scores) to each target site. We then obtained an out-of-frame
    score at a given site by dividing the sum of pattern scores assigned to
    frameshifting deletions by the microhomology score."
    """
    assert(len(seq)>60 and len(seq)<=80)
    seq = seq.upper()
    length_weight=20.0
    right=len(seq)-int(left)

    duplRows = []
    for k in reversed(range(2,left)):
        for j in range(left,left+right-k+1): 
            for i in range(0,left-k+1):
                if seq[i:i+k]==seq[j:j+k]:
                    length = j-i
                    dupSeq = seq[i:i+k]
                    duplRows.append( (dupSeq, i, i+k, j, j+k, length) )

    if len(duplRows)==0:
        return 0, 0

    ### After searching out all microhomology patterns, duplication should be removed!! 
    sum_score_3=0
    sum_score_not_3=0

    for i in range(len(duplRows)):
        n=0
        scrap, left_start, left_end, right_start, right_end, length = duplRows[i]

        for j in range(i):
            _, left_start_ref, left_end_ref, right_start_ref, right_end_ref, _ = duplRows[j]

            if (left_start >= left_start_ref) and \
               (left_end <= left_end_ref) and \
               (right_start >= right_start_ref) and \
               (right_end <= right_end_ref) and \
               (left_start - left_start_ref) == (right_start - right_start_ref) and \
               (left_end - left_end_ref) == (right_end - right_end_ref):
                    n+=1

        if n != 0:
            continue

        length_factor = round(1/math.exp(length/length_weight),3)
        num_GC=scrap.count("G")+scrap.count("C")
        score = 100*length_factor*((len(scrap)-num_GC)+(num_GC*2))

        if (length % 3)==0:
            sum_score_3+=score
        elif (length % 3)!=0:
            sum_score_not_3+=score

        mhScore = sum_score_3+sum_score_not_3
        oofScore = ((sum_score_not_3)*100) / (sum_score_3+sum_score_not_3)
    return int(mhScore), int(oofScore)

def parseEncode(dirName):
    seqSegs = defaultdict(dict)
    for fname in glob.glob(dirName+"/*.bed"):
        cell = basename(fname).replace("wgEncodeAwgSegmentationCombined", "").split(".")[0]
        for line in open(fname):
            seq, seg = line.strip().split()
            seqSegs[seq].setdefault(cell, []).append(seg)
    return seqSegs

    
def sumReads(gSeq):
    " return dict with name -> total number of reads obtained "
    ret = {}
    for name, beds in gSeq.iteritems():
        total = 0
        for bed in beds:
            score = bed[1]
            total+= score
        ret[name] = total
    return ret
    
def main():
    headers = ["offtargetSeq________", "readCount", "copyCount", "guideName", "readShare", "cutType", "hitScore", "guideSeq", "gcGuide", "flank100Gc", "ggMotifCount", "isMainPam", "isAltPam", "mhScore", "effScore", "repCount", "chromatinType"]
    print "\t".join(headers)

    seqToSeg = parseEncode("chromatin/annot/")
    #seqToFaireSeq = parseEncode("chromatin/faireSeq/")
    offtargets = parseBeds("guideSeq")
    totals = sumReads(offtargets)
    flanks = parseFlanks("guideSeq/flankSeq/")
    guideSeqs, guideExtSeqs = parseGuides()

    rowsByType = defaultdict(list)
    for name, guideSeq in guideSeqs.iteritems():
        totalReads = totals[name]
        for off in offtargets[name]:
            offSeq, readCount, seqCount, guideName = off
            flankSeq = flanks[offSeq]
            repCount = flankSeq.count("a")+flankSeq.count("c")+flankSeq.count("t")+flankSeq.count("g")
            offSeq = offSeq.upper()
            flankSeq = flankSeq.upper()

            flankGc = gcCont(flankSeq)
            ggCount = flankSeq[60:80].count("GG")
            hitScore = calcHitScore(offSeq[:20], guideSeq[:20])
            gcGuide = gcCont(guideSeq)
            if gcGuide > 70:
                continue

            seq80 = flankSeq[10:90]
            mhScore, oofScore = calcMicroHomolScore(seq80, 47) # breakpoint at -3 before NGG site (at 50)

            seq30 = flankSeq[36:66]
            effScore = int(100*calcDoenchScore(seq30))

            isMainPam = offSeq.endswith("GG")
            isAltPam = offSeq.endswith("GA") or offSeq.endswith("AG")

            #isFaire = seqToFaireSeq.get(offSeq, False)

            # encode annotations
            encCells = seqToSeg[offSeq]
            #encCellStrs = ["%s=%s" % (cell, ",".join(types)) for cell, types in encCells.items()]
            #encCellStr = "|".join(encCellStrs)
            counter = Counter()
            for cell, types in encCells.iteritems():
                if len(types)==1:
                    counter[types[0]]+=1
            if len(counter)==0:
                encType = "NoAnnot"
            else:
                encType, encCount = counter.most_common(1)[0]

            readShare = float(readCount) / totalReads
            if readShare > 0.01:
                guideType = "strong"
            else:
                guideType = "weak"
            readShare = "%.3f" % (readShare)

            row = list(off)
            row.extend([readShare, guideType, hitScore, flankSeq[40:60], gcGuide, flankGc, ggCount, int(isMainPam), int(isAltPam), mhScore, effScore, repCount, encType])
            rowsByType[guideType].append(tuple(row))

    maxCount = len(rowsByType["strong"])

    for rowsByType, rows in rowsByType.iteritems():
        random.shuffle(rows)
        rows = rows[:maxCount]
        rows.sort(key=operator.itemgetter(1), reverse=True)

        for row in rows:
            row = [str(x) for x in row]
            print "\t".join(row)

#main()
def compSeq(str1, str2):
    " return a string that marks mismatches between str1 and str2 with * "
    s = []
    for x, y in zip(str1, str2):
        if x==y:
            s.append(".")
        else:
            s.append("*")
    return "".join(s)

