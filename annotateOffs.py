# annotate guideseq offtargets with all possible scores we have
import glob, copy, sys, math, operator, random, re, collections, tempfile, subprocess, logging, os, types
import pickle
from collections import defaultdict, Counter
from os.path import basename, join, splitext, isfile, dirname
import glob

import svmlight # install with 'sudo pip install svmlight'

logging.basicConfig(loglevel=logging.INFO)

# for the chari code
import time, gzip, platform

# assignment of activity datasets to genomes
datasetToGenome = {
    "housden2015": "dm3",
    "xu2015Train": "hg19",
    "eschstruth" : "danRer10",
    "doench2014-Hs": "hg19",
    "doench2014-Mm": "mm9",
    "doench2014-CD33Exon2": "hg19",
    "varshney2015": "danRer10",
    "gagnon2014": "danRer10",
    "xu2015": "hg19",
    "ren2015": "dm3",
    "farboud2015": "ce6",
    "schoenig": "hg19",
    "schoenigHs": "hg19",
    "schoenigMm": "mm9",
    "schoenigRn": "rn5",
    "schoenigMm-part2": "mm9",
    "concordet2-Hs": "hg19",
    "concordet2-Mm": "mm9",
    "concordet2-Rn": "rn5",
    "concordet2": "hg19",
    "morenoMateos2015": "danRer7",
    "morenoMateos2015": "danRer7",
    "alenaNonYuvia" : "danRer10",
    "alenaOthers" : "danRer10",
    "alenaPerrine" : "danRer10",
    "alenaYuvia" : "danRer10",
    "chari2015Train": "hg19"
}

# support scoring types when reading .scores.tab files
scoreTypes = ["wang", "wangOrig", "doench", "ssc", "chariRank", "chariRaw", "crisprScan", 'drsc', "fusi"]

def getScoreTypes():
    return scoreTypes

def iterTsvRows(inFile, fieldSep="\t", isGzip=False, skipLines=None, \
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

def calcMitGuideScore(hitSum):
    """ Sguide defined on http://crispr.mit.edu/about 
    Input is the sum of all off-target hit scores. Returns the specificity of the guide.
    """
    score = 100 / (100+hitSum)
    score = int(round(score*100))
    return score

def calcMitGuideScore_offs(guideSeq, otSeqs, maxMm=None, minHitScore=None, minAltHitScore=None):
    " calc mit spec score given a guide and a list of off-target sequences "
    if maxMm:
        newOtSeqs = []
        for otSeq in otSeqs:
            mmCount, diffLogo = countMms(otSeq, guideSeq)
            if mmCount <= maxMm:
                newOtSeqs.append(otSeq)
        otSeqs = newOtSeqs
            
    # split seqs into main and alt
    altOts = []
    mainOts = []
    for ot in otSeqs:
        assert(len(ot)==23)
        if ot.endswith("AG") or ot.endswith("GA"):
            altOts.append(ot)
        else:
            assert(ot.endswith("GG"))
            mainOts.append(ot)

    # calc and filter hit scores
    mainHitScores = [calcHitScore(guideSeq, ot) for ot in mainOts]
    mainHitScores = [h for h in mainHitScores if h > minHitScore]
    altHitScores  = [calcHitScore(guideSeq, ot) for ot in altOts]
    altHitScores  = [h for h in altHitScores if h > minAltHitScore]

    mainHitScores.extend(altHitScores)

    scoreSum = sum(mainHitScores)
    return calcMitGuideScore(scoreSum)

def calcHitScore(string1,string2, startPos=0):
    """ see 'Scores of single hits' on http://crispr.mit.edu/about 
    startPos can be used to feed sequences longer than 20bp into this function
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

def parseScores(fname):
    " parse a key-val tab-sep file, val is float, return as list (key, val) "
    data = []
    for line in open (fname):
        fs = line.strip().split()
        seq, score = fs
        seq = seq.strip('"')
        score = score.strip('"')
        score = float(score)
        data.append((seq, float(score)))
    return data

def parseSvmOut(fname):
    " parse R SVM output file, return as dict seq -> score "
    data = {}
    for line in open (fname):
        fs = line.strip().split()
        seq, score = fs
        seq = seq.strip('"')
        score = score.strip('"')
        data[seq] = float(score)
    return data

def writeDict(d, fname):
    " write dict as a tab file "
    ofh = open(fname, "w")
    for k, v in d.iteritems():
        if type(v)==types.TupleType:
            ofh.write("%s\t%s\n" % (k, "\t".join([str(x) for x in v])))
        else:
            ofh.write("%s\t%s\n" % (k, str(v)))
    ofh.close()

def readDictList(fname, isFloat=False):
    " read tab-sep file into a defaultdict(list) "
    if not isfile(fname):
        logging.warn("%s does not exist. Returning empty dict" % fname)
        return {}

    logging.info("Reading %s" %fname)
    data = defaultdict(list)
    for line in open(fname):
        if line.startswith("#"):
            continue
        fs = line.rstrip("\n").split("\t")
        if len(fs)==2:
            k, v = fs
            if isFloat:
                v = float(v)
        else:
            k = fs[0]
            v = tuple(fs[1:])
            if isFloat:
                v = tuple([float(x) for x in v])
        data[k].append(v)
    return data

def readDict(fname, isFloat=False):
    " read dict from a tab sep file "
    if not isfile(fname):
        logging.warn("%s does not exist. Returning empty dict" % fname)
        return {}

    logging.info("Reading %s" %fname)
    data = {}
    for line in open(fname):
        fs = line.rstrip("\n").split("\t")
        if len(fs)==2:
            k, v = fs
            if isFloat:
                v = float(v)
        else:
            k = fs[0]
            v = tuple(fs[1:])
            if isFloat:
                v = tuple([float(x) for x in v])
        data[k] = v
    return data

def calcSvmEffScores(seqs):
    """
    returns the SVM-calculated efficiency scores from the Wang/Sabatini/Lander paper
    """
    writeSvmRows(seqs, "/tmp/temp.txt")
    cmd = "cd wangSabatiniSvm/; R --slave --no-save -f scorer.R --args /tmp/temp.txt /tmp/temp.out"
    assert(os.system(cmd)==0)
    return parseSvmOut("/tmp/temp.out")

svmScores = None

def lookupSvmScore(seq):
    " retrieve svm scores from svmScores.tab "
    assert(len(seq)==20)
    seq = seq.upper()
    global svmScores
    if svmScores==None:
        svmScores = readDict("svmScores.tab", isFloat=True)
    return svmScores[seq]

def writeSvmRows(seqs, fname):
    """ write the seqs in wang/sabatini SVM format to a file
    #>>> writeSvmRows(["ATAGACCTACCTTGTTGAAG"])
    """
    tmpFile = open(fname, "w")
    #tmpFile = tempfile.NamedTemporaryFile(prefix="svmR")
    for row in iterSvmRows(seqs):
        tmpFile.write("\t".join([str(x) for x in row]))
        tmpFile.write("\n")
    tmpFile.close()

def iterSvmRows(seqs):
    """ calculate the SVM score from the Wang/Sabatini/Lander paper 
    >>> list(iterSvmRows(["ATAGACCTACCTTGTTGAAG"]))
    [['SEQ', 'BP1A', 'BP1C', 'BP1T', 'BP1G', 'BP2A', 'BP2C', 'BP2T', 'BP2G', 'BP3A', 'BP3C', 'BP3T', 'BP3G', 'BP4A', 'BP4C', 'BP4T', 'BP4G', 'BP5A', 'BP5C', 'BP5T', 'BP5G', 'BP6A', 'BP6C', 'BP6T', 'BP6G', 'BP7A', 'BP7C', 'BP7T', 'BP7G', 'BP8A', 'BP8C', 'BP8T', 'BP8G', 'BP9A', 'BP9C', 'BP9T', 'BP9G', 'BP10A', 'BP10C', 'BP10T', 'BP10G', 'BP11A', 'BP11C', 'BP11T', 'BP11G', 'BP12A', 'BP12C', 'BP12T', 'BP12G', 'BP13A', 'BP13C', 'BP13T', 'BP13G', 'BP14A', 'BP14C', 'BP14T', 'BP14G', 'BP15A', 'BP15C', 'BP15T', 'BP15G', 'BP16A', 'BP16C', 'BP16T', 'BP16G', 'BP17A', 'BP17C', 'BP17T', 'BP17G', 'BP18A', 'BP18C', 'BP18T', 'BP18G', 'BP19A', 'BP19C', 'BP19T', 'BP19G', 'BP20A', 'BP20C', 'BP20T', 'BP20G'], ['ATAGACCTACCTTGTTGAAG', 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1]]
    """
    offsets = {"A":0,"C":1,"T":2,"G":3}
    # construct and write header
    headers = ["SEQ"]
    fields = []
    for i in range(1, 21):
        for n in ["A", "C", "T", "G"]:
            fields.append("BP"+str(i)+n)
    headers.extend(fields)
    yield headers

    for seq in seqs:
        row = []
        row.extend([0]*80)
        for pos, nucl in enumerate(seq):
            nuclOffset = offsets[nucl]
            row[pos*4+nuclOffset] = 1
        assert(len(seq)==20)
        row.insert(0, seq)
        yield row

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
    " parse a fasta file, reture dict id -> seq "
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

def parseFastaAsList(fileObj):
    " parse a fasta file, return list (id, seq) "
    seqs = []
    parts = []
    seqId = None
    for line in fileObj:
        line = line.rstrip("\n")
        if line.startswith(">"):
            if seqId!=None:
                seqs.append( (seqId, "".join(parts)) )
            seqId = line.lstrip(">")
            parts = []
        else:
            parts.append(line)
    if len(parts)!=0:
        seqs.append( (seqId, "".join(parts)) )
    return seqs

def iterFastaSeqs(fileObj):
    " parse a fasta file, yield (id, seq) "
    parts = []
    seqId = None
    for line in fileObj:
        line = line.rstrip("\n")
        if line.startswith(">"):
            if seqId!=None:
                yield (seqId, "".join(parts))
            seqId = line.lstrip(">")
            parts = []
        else:
            parts.append(line)
    if len(parts)!=0:
        yield (seqId, "".join(parts))

def parseFlanks(faDir):
    seqToFlank = {}
    for fname in glob.glob(join(faDir, "*.fa")):
        for shortSeq, seq in parseFasta(open(fname)).iteritems():
            assert(shortSeq not in seqToFlank)
            assert(len(seq)==100)
            seqToFlank[shortSeq] = seq
    return seqToFlank

def gcCont(seq):
    assert(len(seq)==20)
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

binDir = "../crispor/bin/Darwin"
baseDir = "../crispor/"

def calcSscScores(seqs):
    """ calc the SSC scores from the paper Xu Xiao Chen Li Meyer Brown Lui Gen Res 2015 
    >>> calcSscScores(["AGCAGGATAGTCCTTCCGAGTGGAGGGAGG"])
    {'AGCAGGATAGTCCTTCCGAGTGGAGGGAGG': 0.182006}
    """
    assert(len(seqs)!=0) # need at least one sequence
    strList = []
    for s in seqs:
        assert(len(s)==30)
        strList.append("%s 0 0 + dummy" % s)
    sscIn = "\n".join(strList)

    # ../../Darwin/SSC -i /dev/stdin  -o /dev/stdout -l 30 -m matrix/human_mouse_CRISPR_KO_30bp.matrix 
    # AGCAGGATAGTCCTTCCGAGTGGAGGGAGG  187 216 -   MYC_exon3_hg19
    # AGCAGGATAGTCCTTCCGAGTGGAGGGAGG  0 0 -   t
    # AGCAGGATAGTCCTTCCGAGTGGAGGGAGG  187 216 -   MYC_exon3_hg19  0.182006
    sscPath = join(binDir, "SSC")
    matPath = join(baseDir, "bin", "src", "SSC0.1", "matrix", "human_mouse_CRISPR_KO_30bp.matrix")
    cmd = [sscPath, "-i", "/dev/stdin", "-o", "/dev/stdout", "-l", "30", "-m", matPath]
    stdout, stderr = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE).communicate(sscIn)
    scores = {}
    i = 0
    for lineIdx, line in enumerate(stdout.split("\n")):
        fs = line.split()
        if "Processing failed" in line:
            raise Exception("SSC returned error, line %d" % lineIdx)
        seq, score = fs[0], float(fs[-1])
        scores[seq] = score
        lineIdx += 1
        if lineIdx==len(seqs):
            break
    return scores

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
    
def outputOldTable():
    " not used anymore, was a try to use ML to predict OT count "
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

def compSeq(str1, str2):
    " return a string that marks mismatches between str1 and str2 with * "
    s = []
    for x, y in zip(str1, str2):
        if x==y:
            s.append(".")
        else:
            s.append("*")
    return "".join(s)

def removeOneNucl(seq):
    " construct 20 sequences, each one with one nucleotide deleted "
    for i in range(0, 20):
        yield seq[:i]+seq[i+1:]

def countMms(string1, string2):
    """ count mismatches between two strings, return mmCount, diffLogo 
    uses only the first 20 nucleotides
    >>> countMms("CCTGCCTCCGCTCTACTCACTGG", "TCTGCCTCCTTTATACTCACAGG")
    (4, '*........**.*.......')
    """
    if len(string1)==23:
        string1 = string1[:20]
        string2 = string2[:20]

    mmCount = 0
    string1 = string1.upper()
    string2 = string2.upper()
    diffLogo = []
    for pos in range(0, len(string1)):
        if string1[pos]!=string2[pos]:
            mmCount+=1
            diffLogo.append("*")
        else:
            diffLogo.append(".")

    return mmCount, "".join(diffLogo)

def findGappedSeqs(guideSeq, offtargetSeq, minMm=9999):
    """ return list of gapped versions of otSeq with lower mismatch count than mmCount
    the mismatch count of the gapped sequences have to be at least lower than
    minMm otherwise nothing will be returned.0
    >>> findGappedSeqs("AATAGC", "AATTAG")
    (1, ['AATAG(C)'], ['ATTAG'], ['.*...'])
    >>> findGappedSeqs("AATAGC", "AATTAG")
    (1, ['AATAG(C)'], ['ATTAG'], ['.*...'])
    >>> findGappedSeqs("CGTACA", "AAGTCA")
    (1, ['CGT(A)CA'], ['AGTCA'], ['*....'])
    >>> findGappedSeqs("TGGATGGAGGAATGAGGAGT", "GAGGATGGGGAATGAGGAGT")
    (1, ['TGGATGG(A)GGAATGAGGAGT'], ['AGGATGGGGAATGAGGAGT'], ['*..................'])
    >>> findGappedSeqs("TGGATGGAGGAATGAGGAGT", "GAGGATGGGGAATGAGGAGT")
    (1, ['TGGATGG(A)GGAATGAGGAGT'], ['AGGATGGGGAATGAGGAGT'], ['*..................'])
    >>> findGappedSeqs("AATAGC", "AATTAG", 1)
    (1, [], [], [])
    """
    # AATAGC
    # AATTAG
    # becomes:
    # AATAG
    # AATAG 
    # and:
    # CGTACA 
    # AAGTCA -> 4 mismatches
    #     **
    # best gapped version is 
    # CGTAC 
    # AGTCA -> 1 mismatch
    # remove either first or last bp from guide
    offtarget1 = offtargetSeq[1:]
    offtarget2 = offtargetSeq[:-1]
    # remove one bp in turn from guideSeq
    guideSeqs = defaultdict(list)
    logos = defaultdict(list)
    gapPos = defaultdict(list)
    otSeqs = defaultdict(list)
    for i, gappedGuideSeq in enumerate(removeOneNucl(guideSeq)):
        for offtargetSeq in [offtarget1, offtarget2]:
            mm, diffLogo = countMms(offtargetSeq, gappedGuideSeq)
            if mm < minMm:
                minMm = mm
                gapPos[mm].append(i)
                logos[mm].append(diffLogo)
                guideSeqs[mm].append(guideSeq[:i]+"("+guideSeq[i]+")"+guideSeq[i+1:])
                otSeqs[mm].append(offtargetSeq)
    return minMm, guideSeqs[minMm], otSeqs[minMm], logos[minMm]

def countMmsAndLogo(string1, string2):
    " count mismatches between two strings also return mismatch ASCII logo"
    mmCount = 0
    string1 = string1.upper()
    string2 = string2.upper()
    diffLogo = []
    for pos in range(0, len(string1)):
        if string1[pos]!=string2[pos]:
            mmCount+=1
            diffLogo.append("*")
        else:
            diffLogo.append(".")

    return mmCount, "".join(diffLogo)

def parseOfftargets(fname, maxMismatches, onlyAlt, validPams):
    """ parse the annotated validated off-target table and return as dict
    guideSeq -> otSeq -> modifFreq and another dict guideName -> guideSeq
    """
    otScores = defaultdict(dict)
    guideSeqs = dict()
    print "parsing %s" % fname
    skipCount = 0
    for row in iterTsvRows(fname):
        #print fname, int(row.mismatches), maxMismatches
        if int(row.mismatches)>maxMismatches:
            #print "skip", row
            skipCount += 1
            continue
        if validPams!=None and not row.otSeq[-2:] in validPams:
            print "not using off-target %s/%s, PAM is not NGG/NGA/NAG" % (row.name, row.otSeq)
            continue

        guideSeqs[row.name] = row.guideSeq
        if onlyAlt and not row.otSeq[-2:] in ["AG", "GA"]:
            continue
        otScores[row.guideSeq][row.otSeq] = float(row.readFraction)
    print "Skipped %d rows with more than %d mismatches" % (skipCount, maxMismatches)
    return otScores, guideSeqs

def parseRawOfftargets(inFname, removeCellLine=True):
    """ parse the raw list of off-targets, in the format of offtargets.tsv.
    returns list of rows and a dict guideName -> guideSeq 
    """
    targetSeqs = {}
    inRows = []
    for row in iterTsvRows(inFname):
        if removeCellLine:
            # by removing the prefix before /, treat Kim's two cell lines as one experiment
            study = row.name.split("_")[0].split("/")[0]
        if row.type=="on-target":
            targetSeqs[row.name] = row.seq
        else:
            inRows.append(row)
    return inRows, targetSeqs

def parseOfftargetsWithNames(fname, maxMismatches, onlyAlt, validPams, useOtNames=False):
    """ parse the annotated validated off-target table and return as dict
    (guideSeq, guideName) -> otSeq -> modifFreq and another dict guideName -> guideSeq
    """
    otScores = defaultdict(dict)
    guideSeqs = dict()
    print "parsing %s" % fname
    skipCount = 0
    for row in iterTsvRows(fname):
        #print fname, int(row.mismatches), maxMismatches
        if int(row.mismatches)>maxMismatches:
            #print "skip", row
            skipCount += 1
            continue
        if validPams!=None and not row.otSeq[-2:] in validPams:
            print "not using off-target %s/%s, PAM is not NGG/NGA/NAG" % (row.name, row.otSeq)
            continue

        if onlyAlt and not row.otSeq[-2:] in ["AG", "GA"]:
            continue
        name = row.name
        if useOtNames:
            name = row.otName
        otScores[(name, row.guideSeq)][row.otSeq] = float(row.readFraction)
    print "Skipped %d rows with more than %d mismatches" % (skipCount, maxMismatches)
    return otScores

def parseMit(dirName, guideSeqs):
    " parse the MIT csv files, return a dict with guideSeq -> otSeq -> otScore "
    #fnames= glob.glob(dirName+"/*.csv")
    #print targetSeqs
    data = defaultdict(dict)
    for guideName in guideSeqs:
    #"for fname in fnames:
        guideNameNoCell = guideName.replace("/K562", "").replace("/Hap1","")
        fname = join(dirName, guideNameNoCell+".csv")
        study = guideName.split("_")[0]
        #if study in ignoreStudies:
            #continue
        #if guideName not in targetSeqs:
            #print "MIT off-target data without bench data: %s" % guideName
            #continue
        print "parsing %s" % fname
        if not isfile(fname):
            logging.error("MISSING: %s" % fname)
            continue
        for line in open(fname):
            if line.startswith("guide"):
                continue
            fs = line.split(", ")
            otSeq = fs[4]
            score = fs[6]
            if fs[7]=="True":
                # ontarget
                continue
            guideSeq = guideSeqs[guideName]
            data[guideSeq][otSeq]=float(score)
    return data

def parseCrispor(dirName, guideNames, maxMismatches):
    """ parse crispor output files, return as dict guideSeq -> ot seq -> otScore 
    Also return a dict with guideName -> guideSeq
    guideNames can be None, in which case all data is read.
    """
    print("Parsing CRISPR results from dir %s" % dirName)
    predScores = defaultdict(dict)
    #targetSeqs = {}
    #for fname in glob.glob(dirName+"/*.tsv"):
    #print "XX", guideNames
    if guideNames is None:
        fnames = glob.glob(join(dirName, "*.tsv"))
        guideNames = [splitext(basename(fname))[0] for fname in fnames]

    for guideName in guideNames:
        # remove cell lines for KIM et al
        guideName = guideName.replace("/K562", "").replace("/Hap1","")
        # fix a few typos
        fname = join(dirName, guideName+".tsv")
        print "parsing %s" % fname
        for row in iterTsvRows(fname):
            if int(row.mismatchCount)>maxMismatches:
                #print "too many mismatches", row, maxMismatches
                continue
            if "hap" in row.chrom or "random" in row.chrom or "chrUn" in row.chrom:
                #print "strange chrom, skipping %s" % str(row)
                continue
            guideName = splitext(basename(fname))[0]
            predScores[row.guideSeq][row.offtargetSeq] = float(row.offtargetScore)
            if row.offtargetSeq=="GAATCCTAAATACTCTCCTTCGG":
                print "XX", row.guideSeq, row.offtargetSeq
            #targetSeqs[guideName] = row.guideSeq
    return predScores

compTable = { "a":"t", "A":"T", "t" :"a", "T":"A", "c":"g", "C":"G", "g":"c", "G":"C", "N":"N", "n":"n", 
        "Y":"R", "R" : "Y", "M" : "K", "K" : "M", "W":"W", "S":"S",
        "H":"D", "B":"V", "V":"B", "D":"H", "y":"r", "r":"y","m":"k",
        "k":"m","w":"w","s":"s","h":"d","b":"v","d":"h","v":"b","y":"r","r":"y" }

def revComp(seq):
    newseq = []
    for nucl in reversed(seq):
       newseq += compTable[nucl]
    return "".join(newseq)

def useRanks(vals, doPercent=False, doQuart=False):
    """ replace values in list with their rank 
    >>> useRanks([0.5, 0.1, 1.5, 1.5])
    [1, 0, 2, 2]
    """
    newList = []
    sortedList = list(sorted(vals))
    for x in vals:
        newList.append(sortedList.index(x))

    if doPercent or doQuart:
        list2 = []
        for i, val in enumerate(newList):
            perc = float(val)/len(newList)
            list2.append(perc)
        newList = list2

    if doQuart:
        list2 = []
        for val in newList:
            if val<0.25:
                x = 1
            elif val < 0.5:
                x = 2
            elif val < 0.75:
                x = 3
            elif val < 1.0:
                x = 4
            else:
                assert(False)
            list2.append(x)
        newList = list2

    return newList

def complRna(seq):
    " complement the sequence and translate to RNA "
    newseq = []
    for nucl in seq.upper():
        newseq.append( compTable[nucl].replace("T", "U") )
    return "".join(newseq)

hsuMat = None # dict with (fromNucl, toNucl) -> list of 19 scores
avgFreqs = None # list of 19 scores

def parseHsuMat(fname):
    """ return the hsu 2013 matrix as a dict rnaNucl -> dnaNucl -> list of scores and a list of 19 averages
    #>>> parseHsuMat("./hsu2013/fig2cData.txt")
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

    # normalize
    # "Each frequency was normalized to range from 0 to 1, such that f = (f-fmin) / (fmax-fmin)"
    normAvgs = [(a-min(avgs)) / (max(avgs)-min(avgs)) for a in avgs]
    normMat = {}
    for key, freqs in hsuMat.iteritems():
        normMat[key] = [( f - minMat) / (maxMat - minMat) for f in freqs]
    normAvgs[4] = 0.001
    assert(min(normAvgs)!=0.0)
    return normMat, normAvgs

def calcHsuSuppScore(guideSeq, otSeq, baseDir="./"):
    """ calculate the score described on page 17 of the Hsu et al 2013 supplement PDF
    >>> calcHsuSuppScore("AGTCCGAGCAGAAGAAGAA","AGTCCGAGCAGAAGAAGAG")
    0.4509132855355929
    >>> calcHsuSuppScore("TGTCCGAGCAGAAGAAGAA","AGTCCGAGCAGAAGAAGAA")
    0.007929899123452079
    >>> calcHsuSuppScore("AGTCCGAGCAGAAGAAGAA","AGTCAGAACAGAAGAACAA")
    3.69683024017458e-08
    >>> countMmsAndLogo("AGTCCGAGCAGAAGAAGAA","AGTCAGAACAGAAGAACAA")
    (3, '....*..*........*..')
    """
    assert(len(otSeq)==19)# Hsu ignores pos 0
    assert(len(guideSeq)==19)
    global hsuMat
    global avgFreqs
    if hsuMat is None:
        matFname = baseDir+"hsu2013/fig2cData.txt"
        hsuMat, avgFreqs = parseHsuMat(matFname)

    rnaSeq = complRna(guideSeq)
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
        #print rnaSeq, guideSeq, i, rnaNucl, dnaNucl, f, g
        score *= f * g

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

    return score

def seqsToVecs(seqs, startPos=0, endPos=None):
    vecs = []
    for s in seqs:
        ensPos = len(s)
        vecs.append(seqToVec(s[startPos:endPos]))
    return vecs

def seqToVec(seq, offsets={"A":0,"C":1,"G":2,"T":3}):
    """ convert a 20bp sequence to a 4*20 0/1 vector 
    >>> seqToVec("AAAAATTTTTGGGGGCCCCC")
    [1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0]
    """
    #assert(len(seq)==20)
    row = [0]*len(seq)*4
    for pos, nucl in enumerate(seq):
        nucl = nucl.upper()
        if nucl in offsets:
            nuclOffset = offsets.get(nucl)
        else:
            print seq
            nuclOffset = offsets["other"]
        row[pos*len(offsets)+nuclOffset] = 1
    return row
    
def vecToSeqDicts(coefs):
    " convert a list of 80 floats to 20 dictionaries with A/C/T/G -> float "
    freqs = []
    for i in range(0,20):
        charFreqs = {}
        for nucl, x in zip("ACGT", range(0,4)):
            freq = coefs[i*4+x]
            if freq==0.0:
                continue
            charFreqs[nucl] = freq
        freqs.append(charFreqs)
    return freqs

def printCoef(coefs):
    " print the 80 floats from a seqtovec operation in a nicer way "
    for i in range(0,20):
        row = []
        for nucl, x in zip("ACGT", range(0,4)):
            row.append("%s=%f" % (nucl, coefs[i*4+x]))
        print " ".join(row)

def alnToVec(guideSeq, otSeq):
    """ encode the mismatches/matches between two 20bp sequences into a 3*20 vector 
    the three parts per basepair are:
    1 - is a G->T or A-C mismatch?
    2 - is a T-T, G-G, A-A match or C->T, T-G, C->A, G->A, A->G mismatch?
    3 - is a T->C or C->C match?
    >>> alnToVec("AAGTCCGAGCAGAAGAAGAA","AAGTCCCAGCAGAGGAAGCA")
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0]
    """
    indexPos = {
        "GT" : 0,
        "AC" : 0,
        "TT" : 1,
        "AT" : 1,
        "TA" : 1,
        "CG" : 1,
        "GG" : 1,
        "GC" : 1,
        "CT" : 1,
        "TG" : 1,
        "CA" : 1,
        "GA" : 1,
        "AA" : 1,
        "AG" : 1,
        "TC" : 2,
        "CC" : 2
    }
    #vec = [ 0.0 ] * 60
    vec = [ 0.0 ] * 20
    for i in range(0, 20):
        subIdx = indexPos[ guideSeq[i]+otSeq[i] ]
        #vec[i*3+subIdx] = 1.0
        if guideSeq[i]!=otSeq[i]:
            vec[i] = 1.0
    return vec

def seqsToChariVecs(seqs):
    """ partially copied from generateSVMFile.FASTA.py
    >>> seqsToChariVecs(["CTTCTTCAAGGTAACTGCAGAGG"])
    [(0, [(11.0, 0.0), (12.0, 0.0), (13.0, 1.0), (14.0, 0.0), (21.0, 0.0), (22.0, 1.0), (23.0, 0.0), (24.0, 0.0), (31.0, 0.0), (32.0, 1.0), (33.0, 0.0), (34.0, 0.0), (41.0, 0.0), (42.0, 0.0), (43.0, 1.0), (44.0, 0.0), (51.0, 0.0), (52.0, 1.0), (53.0, 0.0), (54.0, 0.0), (61.0, 0.0), (62.0, 1.0), (63.0, 0.0), (64.0, 0.0), (71.0, 0.0), (72.0, 0.0), (73.0, 1.0), (74.0, 0.0), (81.0, 0.0), (82.0, 0.0), (83.0, 0.0), (84.0, 1.0), (91.0, 0.0), (92.0, 0.0), (93.0, 0.0), (94.0, 1.0), (101.0, 1.0), (102.0, 0.0), (103.0, 0.0), (104.0, 0.0), (111.0, 1.0), (112.0, 0.0), (113.0, 0.0), (114.0, 0.0), (121.0, 0.0), (122.0, 1.0), (123.0, 0.0), (124.0, 0.0), (131.0, 0.0), (132.0, 0.0), (133.0, 0.0), (134.0, 1.0), (141.0, 0.0), (142.0, 0.0), (143.0, 0.0), (144.0, 1.0), (151.0, 0.0), (152.0, 0.0), (153.0, 1.0), (154.0, 0.0), (161.0, 0.0), (162.0, 1.0), (163.0, 0.0), (164.0, 0.0), (171.0, 1.0), (172.0, 0.0), (173.0, 0.0), (174.0, 0.0), (181.0, 0.0), (182.0, 0.0), (183.0, 1.0), (184.0, 0.0), (191.0, 0.0), (192.0, 0.0), (193.0, 0.0), (194.0, 1.0), (201.0, 1.0), (202.0, 0.0), (203.0, 0.0), (204.0, 0.0), (211.0, 0.0), (212.0, 0.0), (213.0, 0.0), (214.0, 1.0)])]
    """
    vecs = []
    for seq in seqs:
        vec = []
        # end index
        for pos in range(0, 21):
            for nuclIdx, char in enumerate("GTCA"):
                vec.append( (float("%d%d" % (pos+1, nuclIdx+1)), float(seq[pos]==char)) )
        vecs.append( (0, vec) )
    return vecs

def calcChariScoresDirect(seqs, baseDir="."):
    """ calc chari scores without an external process using the svmlight module 
    >>> calcChariScoresDirect(["CTTCTTCAAGGTAACTGCAGAGG"])
    [0.5494762125374468]
    """
    chariDir = baseDir+"/sgRNA.Scorer.1.0/"
    modelFname = chariDir+'293T.HiSeq.SP.Nuclease.100.SVM.Model.txt'
    vecs = seqsToChariVecs(seqs)

    model = svmlight.read_model(modelFname)
    predictions = svmlight.classify(model, vecs)
    return predictions

def calcChariScores(seqs, baseDir="."):
    " return dict with chari 2015 scores, returns seq -> (rawScore, relRank), seqs is a list of 23mers "
    # this is mostly copied from scoreMySites.py in the Chari2015 package
    # handles
    #infile = args[0]
    #organism = args[1]
    #species = args[2]
    #header = args[3]
    organism = "Hg"
    species = "SP"
    chariDir = baseDir+"/sgRNA.Scorer.1.0/"
    header = "/tmp/chariScores"

    # know which files to use
    if species=='SP' and organism=='Hg':
            model = chariDir+'293T.HiSeq.SP.Nuclease.100.SVM.Model.txt'
            dist = chariDir+'Hg19.RefFlat.Genes.75bp.NoUTRs.SPSites.SVMOutput.txt'
            pam = 'NGG'
    elif species=='SP' and organism=='Mm':
            model = chariDir+'293T.HiSeq.SP.Nuclease.100.SVM.Model.txt'
            dist = chariDir+'Mm10.RefFlat.Genes.75bp.NoUTRs.SPSites.SVMOutput.txt'
            pam = 'NGG'
    elif species=='ST1' and organism=='Hg':
            model = chariDir+'293T.HiSeq.ST1.Nuclease.100.V2.SVM.Model.txt'
            dist = chariDir+'Hg19.RefFlat.Genes.75bp.NoUTRs.ST1Sites.SVMOutput.txt'
            pam = 'NNAGAAW'

    elif species=='ST1' and organism=='Mm':
            model = chariDir+'293T.HiSeq.ST1.Nuclease.100.V2.SVM.Model.txt'
            dist = chariDir+'Mm10.RefFlat.Genes.75bp.NoUTRs.ST1Sites.SVMOutput.txt'
            pam = 'NNAGAAW'	
    else:
            raise Exception, "Invalid selection! Choose Hg/Mm and SP/ST1"


    # file names
    gRNAFile      = header + '.putativeGRNASequences.fasta'
    svmInputFile  = header + '.SVMInput.txt'
    svmOutputFile = header + '.SVMOutput.txt'
    finalOutput   = header + '.FinalOutput.txt'

    ## first generate the sites
    ##print 'Time: ' + str(time.ctime())
    ##print 'Generating putative gRNA sites...'
    ##runID = 'python identifyPutativegRNASites.py ' + infile + ' ' + pam + ' ' + gRNAFile
    ##p = subprocess.Popen(runID,shell=True)
    ##p.communicate()
    # added by max
    ofh = open(gRNAFile, "w")
    for s in seqs:
        ofh.write(">%s\n%s\n" % (s, s))
    ofh.close()

    # next generate the SVM input file
    print 'Time: ' + str(time.ctime())
    print 'Generating SVM input file from gRNA sequences...'
    runSVMGen = 'python %sgenerateSVMFile.FASTA.py ' % chariDir + gRNAFile + ' ' + svmInputFile
    p = subprocess.Popen(runSVMGen,shell=True)
    p.communicate()

    # run the SVM
    print 'Time: ' + str(time.ctime())
    print 'Running classification using SVM-Light'

    # SVM Classify
    currPlatform = platform.system()
    binPath = dirname(join(__file__))+"/bin/%s/svm_classify" % currPlatform

    runSVMClassify = binPath+' -v 0 ' + svmInputFile + ' ' + model + ' ' + svmOutputFile

    p = subprocess.Popen(runSVMClassify,shell=True)
    p.communicate()

    # write the final outputHeader
    print 'Time: ' + str(time.ctime())
    print 'Converting scores to ranks based on global ' + species + ' score distribution'
    runMakeTable = 'python %smakeFinalTable.py ' % chariDir + gRNAFile + ' ' + svmOutputFile + ' ' + dist + ' ' + finalOutput
    p = subprocess.Popen(runMakeTable,shell=True)
    p.communicate()

    # added by max
    scores = dict()
    for line in open(finalOutput):
        if line.startswith("SeqID"):
            continue
        fs = line.rstrip("\n").split("\t")
        scores[fs[1]] = (float(fs[-2]), float(fs[-1]))
    return scores

chariScores = None

def lookupChariScore(seq):
    " retrieve one chari score from chariScores.tab, return tuple (rawScore, relPercRank) "
    if not len(seq)==23:
        print seq
    assert(len(seq)==23)
    seq = seq.upper()
    global chariScores
    if chariScores==None:
        chariScores = readDict("chariScores.tab", isFloat=True)
    return chariScores[seq]

crisprScanScores = None

def lookupCrisprScan(seq):
    " retrieve one crisprScan score from crisprScanScores.tab, return score "
    assert(len(seq)==23)
    seq = seq.upper()
    global crisprScanScores
    if crisprScanScores==None:
        crisprScanScores = readDict("crisprScanScores.tab", isFloat=True)
    score = crisprScanScores.get(seq)
    if score is None:
        print "No crisprScan score for seq %s. Returning 0" % seq
        return 0
    else:
        return float(score)

myClf = None

def calcEffScores(datasetName, seqs):
    " given list of 34mers, return dict with seq -> scoreName -> score "
    global myClf
    if myClf is None:
        myClf = pickle.load(open("out/svm.pickle"))
        
    sscSeqs = [s[-30:] for s in seqs]
    sscScores = calcSscScores(sscSeqs)

    scores = defaultdict(dict)
    for seq in seqs:
        #print seq, len(seq)
        if (len(seq)!=34):
            logging.error( "Seq is not 34 bp long %s, len = %d" %  (seq, len(seq)))
            assert(False)
        scores[seq]["doench"] = calcDoenchScore(seq[:30])
        scores[seq]["ssc"] = sscScores[seq[-30:]]
        scores[seq]["svm"] = 1.0 - lookupSvmScore(seq[4:24])
        chariRaw, chariRank = lookupChariScore(seq[4:27])
        scores[seq]["chariRaw"] = chariRaw
        scores[seq]["chariRank"] = chariRank
        if not skipOof:
            scores[seq]["oof"]  = lookupOofScore(seq[4:27])
        #scores[seq]["crisprScan"]  = lookupCrisprScan(seq[4:27])

        guideSeq = seq[4:24]
        assert(len(guideSeq)==20)
        scores[seq]["finalGc6"] = countFinalGc(guideSeq, 6)
        #scores[seq]["finalGc2"] = countFinalGc(guideSeq, 2)
        scores[seq]["finalGg"] = (guideSeq[-2:]=="GG")
        vec = seqToVec(guideSeq[10:])
        #scores[seq]["myScore"] = myClf.predict(vec)[0]
    return scores

def countFinalGc(seq, lastCount):
    """ return the  GC count of the 20mer sequence in the last few base pairs
    >>> countFinalGc("ATAGACCTACCTTGTTGAAG", 3)
    1.0
    >>> countFinalGc("AAAAAAAAAAAAAAAAAAAA", 10)
    0.0
    """
    assert(len(seq)==20)
    wgc = 0.0
    for nucl in seq[-lastCount:]:
        if nucl == "G" or nucl=="C":
            wgc += 1.0
    return wgc

def convToRankPerc(vec):
    """ given a list of values, return their rank-percentiles as a list in the same order 
    >>> convToRankPerc([1,1,5,3,3])
    [0.2, 0.2, 0.8, 0.6, 0.6]
    """
    sortVec = list(sorted(vec))
    valToPercRank = dict()
    for i in range(0, len(sortVec)):
        percRank = float(i)/len(sortVec)
        valToPercRank[sortVec[i]] = percRank

    percRankVec = []
    for x in vec:
        percRankVec.append(valToPercRank[x])
    return percRankVec

# sqlUcsc hgcentral -e 'select * from blatServers where db="mm9"'
# rsync -avp hgdownload.soe.ucsc.edu::gbdb/dm3/dm3.2bit ./ --progress --partial
blatServers = {
    "hg19": ("blat1a", "17779"),
    "danRer10" : ("blat1c", "17863"),
    "dm3" : ("blat4d", "17791"),
    "ce6" : ("blat4d", "17841"),
    "rn5" : ("blat1b", "17795"),
    "mm9" : ("blat1c", "17779")
}

def extendSeqs(seqs, db, fiveExt, threeExt):
    """ extend seqs for eff score calculations, seqs is a dict id -> seq
    by fiveExt and threeExt into 5' / 3' direction
    return a dict seqId -> list of (newSeq, genomePositionString)
    """
    ret = defaultdict(list)

    # fill the ret object with everything we have in the cache
    cacheFname ="out/blatCache.tab"
    print "Opening %s" % cacheFname
    blatCache = parseBlatCache(cacheFname)
    for seq in seqs:
        if seq in blatCache:
            ret[seq] = blatCache[seq]

    # return what we have if nothing left to do
    if len(set(seqs)-set(ret))==0:
        return ret
    # remove what we got from the cache from seqs
    for seq in ret:
        seqs.remove(seq)

    # now start the blatting procedure
    # write 23mers to fa file
    inNames = set()
    ofh = open("/tmp/temp.fa", "w")
    seqLen = None
    for seqId, seq in seqs.iteritems():
        ofh.write(">%s\n%s\n" % (seqId, seq))
        if seqLen!=None:
            assert(len(seq)==seqLen)
        else:
            seqLen=len(seq)
        inNames.add(seqId)
    ofh.close()

    print "running BLAT, writing to /tmp/temp.bed"
    blatServer, port = blatServers[db]
    cmd = "gfClient %s.soe.ucsc.edu %s /gbdb/%s /tmp/temp.fa /tmp/temp.psl -minScore=20 -nohead -minIdentity=100 -maxIntron=0 -dots=1 ; pslToBed /tmp/temp.psl /tmp/temp.bed" % (blatServer, port, db)
    os.system(cmd)

    matches = defaultdict(list) # seqId -> list of (chrom, start, end, strand)
    for line in open("/tmp/temp.bed"):
        chrom, start, end, name, score, strand = line.split()[:6]
        if "_hap" in chrom or "random" in chrom or "chrUn" in chrom:
            continue
        ##print int(end)-int(start)
        if (int(end)-int(start))!=seqLen:
            #print "SKIP", line
            continue
        matches[name].append( (chrom, int(start), int(end), strand) )

    notFoundNames = inNames - set(matches)
    logging.warn("These sequences were not found with BLAT: %s" % ",".join(notFoundNames))
    #assert( len(seqs) == len(matches) )

    # write matches to temp.bed file
    # SSC needs extension by +7 bp of the end position
    # Doench needs extension -4 of the start and +3 of the end pos
    print "Creating /tmp/tempExt.bed with extended matches"
    ofh = open("/tmp/tempExt.bed", "w")
    positions = []
    for seqId, matchTuples in matches.iteritems():
        if len(matchTuples)>1:
            logging.error("Multiple matches for %s, will require manual selection" % seqId)
            logging.error("%s" % matchTuples)
        for matchTuple in matchTuples:
            chrom, start, end, strand = matchTuple
            if strand=="+":
                start = start - fiveExt
                end = end + threeExt
            else:
                start = start - threeExt
                end = end + fiveExt
            row = [chrom, str(start), str(end), seqId, "0", strand]
            ofh.write("\t".join(row)+"\n")
            positions.append( "%s:%d-%d:%s" % (chrom, start, end, strand))
    ofh.close()

    cmd = "twoBitToFa /gbdb/%s/%s.2bit -bed=/tmp/tempExt.bed /tmp/tempExt.fa" % (db, db)
    os.system(cmd)
    seqs = parseFastaAsList(open("/tmp/tempExt.fa"))
    assert(len(seqs)==len(positions))

    for seqData, pos in zip(seqs, positions):
        seqId, seq = seqData
        ret[seqId].append( (seq, pos) )

    return ret

def addDoenchAndScs(fname):
    """ given tab file with extSeq column, yield (guide, modFreq, scoreType -> score) 
    """
    seqs = []
    freqs = {}
    for row in iterTsvRows(fname):
        seqs.append(row.extSeq)
        modFreq = row.modFreq
        #if modFreq=="None":
            #modFreq = row.geneAbund
        freqs[row.extSeq] = (row.guide, float(modFreq))
    effScores = calcEffScores(seqs)

    #for row in iterTsvRows(fname):
        #yield row.guide, row.modFreq, effScores[row.seq]
    return effScores, freqs

def extend23Mers(seqs, db):
    """ extend 23mers to 34mers for eff score calculations, seqs is a dict id -> seq 
    return a dict seqId -> list of (seq, genomePositionString)
    """
    print "extending %d 23mer sequences on genome" % len(seqs)
    ret = defaultdict(list)

    # fill the ret object with everything we have in the cache
    cacheFname ="out/blatCache.tab"
    print "Opening %s" % cacheFname
    blatCache = parseBlatCache(cacheFname)
    for seqId, seq in seqs.iteritems():
        if seq in blatCache:
            ret[seqId] = blatCache[seq]

    # remove what we got from seqs
    for seqId in ret:
        del seqs[seqId]

    if len(seqs)==0:
        return ret

    print "%d left for BLAT" % len(seqs)
    print seqs

    # write 23mers to fa file
    inNames = set()
    ofh = open("/tmp/temp.fa", "w")
    for seqId, seq in seqs.iteritems():
        ofh.write(">%s\n%s\n" % (seqId, seq))
        inNames.add(seqId)
    ofh.close()

    print "running BLAT, writing to /tmp/temp.bed"
    blatServer, port = blatServers[db]
    cmd = "gfClient %s.soe.ucsc.edu %s /gbdb/%s /tmp/temp.fa /tmp/temp.psl -minScore=20 -nohead -minIdentity=100 -maxIntron=0 -dots=1 ; pslToBed /tmp/temp.psl /tmp/temp.bed" % (blatServer, port, db)
    os.system(cmd)

    matches = defaultdict(list) # seqId -> list of (chrom, start, end, strand)
    for line in open("/tmp/temp.bed"):
        chrom, start, end, name, score, strand = line.split()[:6]
        if "_hap" in chrom or "random" in chrom:
            continue
        ##print int(end)-int(start)
        if (int(end)-int(start))!=23:
            #print "SKIP", line
            continue
        matches[name].append( (chrom, int(start), int(end), strand) )

    notFoundNames = inNames - set(matches)
    logging.warn("These sequences were not found with BLAT: %s" % ",".join(notFoundNames))
    #assert( len(seqs) == len(matches) )

    # write matches to temp.bed file
    # SSC needs extension by +7 bp of the end position
    # Doench needs extension -4 of the start and +3 of the end pos
    print "Creating /tmp/tempExt.bed with extended matches"
    ofh = open("/tmp/tempExt.bed", "w")
    positions = []
    for seqId, matchTuples in matches.iteritems():
        if len(matchTuples)>1:
            logging.error("Multiple matches for %s, will require manual selection" % seqId)
            logging.error("%s" % matchTuples)
        for matchTuple in matchTuples:
            chrom, start, end, strand = matchTuple
            if strand=="+":
                start = start - 4
                end = end + 7
            else:
                start = start - 7
                end = end + 4
            row = [chrom, str(start), str(end), seqId, "0", strand]
            ofh.write("\t".join(row)+"\n")
            positions.append( "%s:%d-%d:%s" % (chrom, start, end, strand))
    ofh.close()

    cmd = "twoBitToFa /gbdb/%s/%s.2bit -bed=/tmp/tempExt.bed /tmp/tempExt.fa" % (db, db)
    os.system(cmd)
    seqs = parseFastaAsList(open("/tmp/tempExt.fa"))
    assert(len(seqs)==len(positions))

    for seqData, pos in zip(seqs, positions):
        seqId, seq = seqData
        ret[seqId].append( (seq, pos) )

    return ret

def extendTabAddContext(fname, db):
    """
    add a column to a tab file with a seq column: new column is called seqExt with 34 mers
    additional column is pos, with chrom:start-end:strand
    Is a NoOp if the ext.tab file already exists
    """
    newFname = fname.replace(".tab", ".ext.tab")
    if isfile(newFname):
        logging.info("not recreating %s, already exists. Delete to recreate" % newFname)
        return newFname
    if db==None:
        return newFname

    seqs = dict()
    for row in iterTsvRows(fname):
        seqs[row.guide] = row.seq
    seqPosDict = extend23Mers(seqs, db)

    ofh = open(newFname, "w")
    ofh.write("\t".join(row._fields)+"\t")
    ofh.write("extSeq\tposition\n")

    for row in iterTsvRows(fname):
        guideName = row.guide
        for seqPos in seqPosDict[guideName]:
            seq, pos = seqPos
            ofh.write("\t".join(row)+"\t"+seq+"\t"+pos+"\n")
    ofh.close()
    print "Wrote result to %s" % ofh.name
    return newFname

def parseSeqScores(datasetName):
    inFname = join("effData/"+datasetName+".scores.tab")
    seqs, scores = [], []
    for row in iterTsvRows(inFname):
        seqs.append(row.seq)
        scores.append(float(row.modFreq))
    return seqs, scores

def parseEffScores(datasetName):
    """ parse an efficiency dataset from the effData/ directory and return 
    a dict seq -> scoreType -> score and dict seq -> (guideName, modFreq)
    >>> parseEffScores("xu2015")
    """
    inFname = join("effData/"+datasetName+".scores.tab")
    print "reading %s" % inFname
    #extFname = extendTabAddContext(inFname, db)
    #scores, freqs = addDoenchAndScs(extFname)

    scores = {}
    freqs  = {}
    for row in iterTsvRows(inFname):
        seq = row.seq
        scores[seq] = {}
        for st in scoreTypes:
            scores[seq][st] = float(row._asdict()[st])

        #scores[seq]["wang"] = scores[seq]["wang"]

        # stay compatible with new files where all scores all 0-100
        if scores[seq]["doench"] > 1.0:
            scores[seq]["doench"] = scores[seq]["doench"] / 100
        if scores[seq]["wangOrig"] > 1.0:
            scores[seq]["wangOrig"] = scores[seq]["wangOrig"] / 100
        if scores[seq]["wang"] > 1.0:
            scores[seq]["wang"] = scores[seq]["wang"] / 100

        guideSeq = row.seq[:20]
        scores[seq]["finalGc6"] = int(countFinalGc(guideSeq, 6)>=4)
        scores[seq]["finalGg"] = int(guideSeq[-2:]=="GG")

        freq = float(row.modFreq)
        # inverse sign -> higher = better
        if datasetName.startswith("xu2015Train"):
            freq = -freq

        freqs[seq] = (row.guide, freq)
        
    assert(len(scores)!=0)
    assert(len(scores)==len(freqs))
    return scores, freqs

def parseOofScores(fname):
    " read the file out/oofScores.tab generated by calcOofScores.py and return a dict 23mer -> oofScore "
    ret = {}
    for row in iterTsvRows(fname):
        ret[row.seq] = float(row.oofScore)
    return ret

oofScores = None

def lookupOofScore(seq):
    " given a 23mer, return its OOF score "
    if not len(seq)==23:
        print seq
    assert(len(seq)==23)
    global oofScores
    if oofScores is None:
        oofScores = parseOofScores("out/oofScores.tab")
    score = oofScores.get(seq, None)
    if score==None:
        print "couldn't find OOF score for %s" % seq
        score = None
    return score
    #return oofScores[seq]

def parseBlatCache(fname):
    " return dict 23mer -> list of (34mer, genomePos) "
    data = defaultdict(list)
    for row in iterTsvRows(fname):
        data[row.seq].append((row.extSeq, row.position))
    return data

def parseAllGuides():
    mask = join("effData/"+"*.tab")
    seqs, scores = [], []
    doneSeqs = set()
    for fname in glob.glob(mask):
        for row in iterTsvRows(fname):
            if row.seq in doneSeqs:
                continue
            seqs.append(row.seq)
            doneSeqs.add(row.seq)
            scores.append(float(row.modFreq))
    return seqs, scores

def writeRow(ofh, row):
    " write list to file as tab-sep row "
    row = [str(x) for x in row]
    ofh.write("\t".join(row))
    ofh.write("\n")

if __name__=="__main__":
    import doctest
    doctest.testmod()
