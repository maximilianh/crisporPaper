# this library re-implements the efficiency scoring functions of these articles:

# - WangSvm: Wang et al, Science 2014, PMID 24336569, no website
# - Doench: Doench et al, Nat Biotech 2014, PMID 25184501, http://www.broadinstitute.org/rnai/public/analysis-tools/sgrna-design
# - CrisprScan: Moreno-Mateos, Nat Meth 2015, PMID 26322839, http://crisprscan.org
# - SCS: Xu et al, Gen Res 2015, PMID 26063738, http://crispr.dfci.harvard.edu/SSC/
# - Chari: Chari et al, PMID 26167643 http://crispr.med.harvard.edu/sgRNAScorer

# the input are 100bp sequences that flank the basepair just 5' of the PAM +/-50bp.
# so 50bp 5' of the PAM, and 47bp 3' of the PAM -> 100bp

# this module uses pipes to feed data into some programs
# If you run too many sequences at once, it may hang. Increase the BUFSIZE parameter in this case.

from subprocess import Popen, PIPE, STDOUT, check_output, CalledProcessError, call
import platform, math, tempfile, bisect, sys, os, logging, types
from os.path import dirname, join, basename, isfile, expanduser

import urllib2
import json

# import numpy as np

# global that points to the crispor 'bin' directory with the external executables
# like libsvm and svmlight
global binDir
binDir = None

if binDir is None:
    binDir = join(".", "bin")

BUFSIZE = 10000000

def setBinDir(path):
    global binDir
    binDir = path

def getBinPath(name):
    """
    get the full pathname of a platform-specific binary, in the bin/ directory relative to this directory
    """
    currPlatform = platform.system()
    #myDir = dirname(join(__file__))
    binPath = join(binDir, currPlatform, name)
    if not isfile(binPath):
        raise Exception("Could not find %s" % binPath)
    return binPath

def seqToVec(seq, offsets={"A":0,"C":1,"G":2,"T":3}):
    """ convert a x bp sequence to a 4 * x 0/1 vector
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

paramsCRISPRscan = [
# converted excel table of logistic regression weights with 1-based positions
('AA',18,-0.097377097),
('TT',18,-0.094424075),('TT',13,-0.08618771),('CT',26,-0.084264893),('GC',25,-0.073453609),
('T',21,-0.068730497),('TG',23,-0.066388075),('AG',23,-0.054338456),('G',30,-0.046315914),
('A',4,-0.042153521),('AG',34,-0.041935908),('GA',34,-0.037797707),('A',18,-0.033820432),
('C',25,-0.031648353),('C',31,-0.030715556),('G',1,-0.029693709),('C',16,-0.021638609),
('A',14,-0.018487229),('A',11,-0.018287292),('T',34,-0.017647692),('AA',10,-0.016905415),
('A',19,-0.015576499),('G',34,-0.014167123),('C',30,-0.013182733),('GA',31,-0.01227989),
('T',24,-0.011996172),('A',15,-0.010595296),('G',4,-0.005448869),('GG',9,-0.00157799),
('T',23,-0.001422243),('C',15,-0.000477727),('C',26,-0.000368973),('T',27,-0.000280845),
('A',31,0.00158975),('GT',18,0.002391744),('C',9,0.002449224),('GA',20,0.009740799),
('A',25,0.010506405),('A',12,0.011633235),('A',32,0.012435231),('T',22,0.013224035),
('C',20,0.015089514),('G',17,0.01549378),('G',18,0.016457816),('T',30,0.017263162),
('A',13,0.017628924),('G',19,0.017916844),('A',27,0.019126815),('G',11,0.020929039),
('TG',3,0.022949996),('GC',3,0.024681785),('G',14,0.025116714),('GG',10,0.026802158),
('G',12,0.027591138),('G',32,0.03071249),('A',22,0.031930909),('G',20,0.033957008),
('C',21,0.034262921),('TT',17,0.03492881),('T',13,0.035445171),('G',26,0.036146649),
('A',24,0.037466478),('C',22,0.03763162),('G',16,0.037970942),('GG',12,0.041883009),
('TG',18,0.045908991),('TG',31,0.048136812),('A',35,0.048596259),('G',15,0.051129717),
('C',24,0.052972314),('TG',15,0.053372822),('GT',11,0.053678436),('GC',9,0.054171402),
('CA',30,0.057759851),('GT',24,0.060952114),('G',13,0.061360905),('CA',24,0.06221937),
('AG',10,0.063717093),('G',10,0.067739182),('C',13,0.069495944),('GT',31,0.07342535),
('GG',13,0.074355848),('C',27,0.079933922),('G',27,0.085151052),('CC',21,0.088919601),
('CC',23,0.095072286),('G',22,0.10114438),('G',24,0.105488325),('GT',23,0.106718563),
('GG',25,0.111559441),('G',9,0.114600681)]

def calcCrisprScanScores(seqs):
    """ input is a 35bp long sequence: 6bp 5', 20bp guide, 3 bp PAM and 6bp 3'
    >>> calcCrisprScanScores(["TCCTCTGGTGGCGCTGCTGGATGGACGGGACTGTA"])
    [0.7738634886289999]
    """
    scores = []
    for seq in seqs:
        assert(len(seq)==35)
        intercept = 0.183930943629
        score = intercept
        for modelSeq, pos, weight in paramsCRISPRscan:
            subSeq = seq[pos-1:pos+len(modelSeq)-1]
            if subSeq==modelSeq:
                score += weight
        scores.append(score)
    return scores

def listToSvml(vec, res):
    """ convert a list of values to a line in svml format line like "+1 1:0.5 2:1.5 ...
    """
    parts = [str(res)]
    for i, val in enumerate(vec):
        parts.append("%d:%d" % (i+1, val))
    return " ".join(parts)

def calcWangSvmScores(seqs):
    """
    Use the wang.model file to score sequences. Input is only the 20bp guide sequence.
    Uses libsvm's svm-predict program, V2.6.
    The results here are off mostly by 1-5% from the results returned by the Wang et al source code.
    I never found out why, there are no parameters for "svm_predict". Should not be due to a version
    difference either, I'm using the same libsvm version as the e1071 R module.
    This is necessary for a web server implementation as e1071 in R cannot read the model from a file.

    The original implementation from the paper can be called with calcWangSvmScoresUsingR()

    See compareWangScores.py:

    The Pearson correlation between both ways to calculate the score is 97%.

    Histogram of the score differences:
    0.000000 ************************************************************ 3074
    0.050000 ********************************* 1674
    0.100000 ************ 612
    0.150000 **** 191
    0.200000 * 52
    0.250000  7
    0.300000  1
    car out/wangDiffs.tsv | cut -f4 | tr -d '-' | grep -v diff | textHistogram stdin stdout -real -binSize=0.05

    >>> calcWangSvmScores(["ATAGACCTACCTTGTTGAAG"])
    [0.395414]
    """
    scores = []
    vecOrder = {"A":0, "C":1, "T":2, "G":3}

    lines = []
    for seq in seqs:
        assert(len(seq)==20)
        vec = seqToVec(seq, offsets=vecOrder)
        lines.append(listToSvml(vec, 0))

    dataIn = "\n".join(lines)
    binPath = getBinPath("svm-predict")
    modelFname = join(binDir, "src", "wangSabatiniSvm", "wang.model")
    cmd = [binPath, "-b", "1", "/dev/stdin", modelFname, "/dev/stdout"]
    proc = Popen(cmd,stdout=PIPE, stdin=PIPE, stderr=None, bufsize=BUFSIZE)
    dataOut = proc.communicate(input=dataIn)[0]

    for line in dataOut.splitlines():
        if line.startswith("labels"):
            continue
        if line.startswith("Accuracy"):
            break
        score = float(line.split()[-1])
        scores.append(score)

    return scores

# DOENCH SCORING 
doenchParams = [
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

def calcDoenchScores(seqs):
    """
    Code reproduced following paper's methods section. Thanks to Daniel McPherson for fixing it.
    Input is a 30mer: 4bp 5', 20bp guide, 3bp PAM, 3bp 5'
    """
    intercept =  0.59763615
    gcHigh    = -0.1665878
    gcLow     = -0.2026259

    scores = []
    for seq in seqs:
        assert(len(seq)==30)
        score = intercept

        guideSeq = seq[4:24]
        gcCount = guideSeq.count("G") + guideSeq.count("C")
        if gcCount <= 10:
            gcWeight = gcLow
        if gcCount > 10:
            gcWeight = gcHigh
        score += abs(10-gcCount)*gcWeight

        for pos, modelSeq, weight in doenchParams:
            subSeq = seq[pos:pos+len(modelSeq)]
            if subSeq==modelSeq:
                score += weight
        expScore = 1.0/(1.0+math.exp(-score))
        scores.append(expScore)

    return scores

def calcSscScores(seqs):
    """ calc the SSC scores from the paper Xu Xiao Chen Li Meyer Brown Lui Gen Res 2015 
    Input is a 30mer, 20bp for the guide, 3bp PAM, 7bp 3' flanking
    >>> calcSscScores(["AGCAGGATAGTCCTTCCGAGTGGAGGGAGG"])
    [0.182006]
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
    sscPath = getBinPath("SSC")
    matPath = join(binDir, "src", "SSC0.1", "matrix", "human_mouse_CRISPR_KO_30bp.matrix")
    cmd = [sscPath, "-i", "/dev/stdin", "-o", "/dev/stdout", "-l", "30", "-m", matPath]
    try:
        stdout, stderr = Popen(cmd, stdin=PIPE, stdout=PIPE, bufsize=BUFSIZE).communicate(sscIn)
    except OSError:
        raise Exception("Cannot run command %s" % " ".join(cmd))
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

    scoreList = []
    # make sure we got a score for each input sequence
    for s in seqs:
        scoreList.append(scores[s])
        
    return scoreList

def seqsToChariSvml(seqs):
    """ partially copied from generateSVMFile.FASTA.py in the Chari et al source code
    >>> seqsToChariSvml(["CTTCTTCAAGGTAACTGCAGA", "CTTCTTCAAGGTAACTGGGGG"])
    '0 13:1 22:1 32:1 43:1 52:1 62:1 73:1 84:1 94:1 101:1 111:1 122:1 134:1 144:1 153:1 162:1 171:1 183:1 194:1 201:1 214:1\\n0 13:1 22:1 32:1 43:1 52:1 62:1 73:1 84:1 94:1 101:1 111:1 122:1 134:1 144:1 153:1 162:1 171:1 181:1 191:1 201:1 211:1'
    """
    vecs = []
    for seq in seqs:
        assert(len(seq)==21)
        vec = []
        # end index
        for pos in range(0, 21):
            for nuclIdx, char in enumerate("GTCA"):
                val = int(seq[pos]==char)
                if val!=0:
                    vec.append( ("%d%d" % (pos+1, nuclIdx+1), val) )
        vecs.append( vec )

    lines = []
    for vec in vecs:
        vec = ["%s:%s" % (x,y) for x,y in vec]
        lines.append("0 "+" ".join(vec))
    return "\n".join(lines)

chariRanges = None

def convertChariToRankPerc(score):
    """
    convert chari score to rank percent using only 2 digits precision. A lot faster than
    the version in the paper.
    """
    global chariRanges
    if chariRanges is None:
        # parse values
        fname = join(binDir, "src", "sgRNA.Scorer.1.0", "Hg19.RefFlat.Genes.75bp.NoUTRs.SPSites.SVMOutput.ranges.txt")
        ranges = open(fname).read().splitlines()
        ranges = [float(x) for x in ranges]

    # use bisection to find the right value
    fastPerc = bisect.bisect(ranges, score)-1

    # the old, slow way
    #fname = join(binDir, "src", "sgRNA.Scorer.1.0", "Hg19.RefFlat.Genes.75bp.NoUTRs.SPSites.SVMOutput.txt")
    #allData = open(fname).read().splitlines()
    #allData = np.array([float(x) for x in allData])
    #slowPerc = 100.0*(allData[allData < score].size / float(allData.size))
    return fastPerc

def calcChariScores(seqs, baseDir="."):
    """ return dict with chari 2015 scores, returns two lists (rawScores, rankPercent)
    input seqs have lengths 21bp: 20 bp guide + 1bp first from PAM
    >>> calcChariScores(["CTTCTTCAAGGTAACTGCAGA", "CTTCTTCAAGGTAACTGGGGG"])
    ([0.54947621, 0.58604487], [80, 81])
    """
    # this is a rewritten version of scoreMySites.py in the Chari2015 suppl files
    chariDir = join(binDir, "src", "sgRNA.Scorer.1.0")
    modelFname = join(chariDir,'293T.HiSeq.SP.Nuclease.100.SVM.Model.txt')
    dataIn = seqsToChariSvml(seqs)
    #print repr(dataIn)

    #tempFh = tempfile.NamedTemporaryFile()
    tempFname = tempfile.mktemp()
    #tempFh = open("temp3.txt", "w")
    tempFh = open(tempFname, "w")
    tempFh.write(dataIn+"\n")
    tempFh.close()
    #tempFname = tempFh.name
    #tempFh.close()

    #outTempFh = tempfile.NamedTemporaryFile()
    #outName = outTempFh.name
    outName = tempfile.mktemp()

    svmlPath = getBinPath("svm_classify")
    cmd = [svmlPath, "-v", "0", tempFname, modelFname, outName]
    try:
        proc = call(cmd)
    except CalledProcessError:
        raise Exception("Could not run command '%s'" % (" ".join(cmd)))

    #print " ".join(cmd)
    dataOut = open(outName).read()
    os.remove(outName)
    os.remove(tempFname)

    scores = []
    ranks = []
    for line in dataOut.splitlines():
        score = float(line)
        scores.append(score)
        ranks.append(convertChariToRankPerc(score))
    return scores, ranks

    #cmd = svmlight.classify(model, vecs)
    #return scores

def writeDict(d, fname):
    " write dict as a tab file "
    ofh = open(fname, "w")
    for k, v in d.iteritems():
        if type(v)==types.TupleType:
            ofh.write("%s\t%s\n" % (k, "\t".join([str(x) for x in v])))
        else:
            ofh.write("%s\t%s\n" % (k, str(v)))
    ofh.close()

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

def sendFusiRequest(seqs):
    """ obtain the fusi score as calculated by Fusi et al's webservice
    >>> sendFusiRequest([ "GGGAGGCTGCTTTACCCGCTGTGGGGGCGC", "GGGAGGCTGCTTTACCCGCTGTGGGGGCGC"])
    [0.599829103280437, 0.599829103280437]
    """
    cacheFname = "out/fusiCache.tab"
    scoreCache = readDict(cacheFname)

    # select sequences that are not in cache yet
    reqSeqs = set()
    for s in seqs:
        if not s in scoreCache:
            reqSeqs.add(s)

    if len(reqSeqs)==0:
        scoreList = [scoreCache[s] for s in seqs]
        return scoreList

    keyFname = expanduser("~/.fusiKey.txt")
    if not isfile(keyFname):
        raise Exception("No ~/.fusiKey.txt file found. Request an API key from azimuth@microsoft.com, write it into this file (single line) and retry")

    api_key = open(keyFname, "r").read().strip()
    paramList = [ (seq, "-1", "-1") for seq in reqSeqs]
                        #"Values": [ [ "GGGAGGCTGCTTTACCCGCTGTGGGGGCGC", "-1", "-1" ] ]
    data =  {

            "Inputs": {

                    "input1":
                    {
                        "ColumnNames": ["sequence", "cutsite", "percentpeptide"],
                        "Values": paramList,
                    },        },
                "GlobalParameters": {
    }
        }

    body = str.encode(json.dumps(data))

    url = 'https://ussouthcentral.services.azureml.net/workspaces/ee5485c1d9814b8d8c647a89db12d4df/services/c24d128abfaf4832abf1e7ef45db4b54/execute?api-version=2.0&details=true'
    headers = {'Content-Type':'application/json', 'Authorization':('Bearer '+ api_key)}

    req = urllib2.Request(url, body, headers)

    try:
        response = urllib2.urlopen(req)

        # If you are using Python 3+, replace urllib2 with urllib.request in the above code:
        # req = urllib.request.Request(url, body, headers) 
        # response = urllib.request.urlopen(req)

        dataList = json.loads(response.read())["Results"]["output2"]["value"]["Values"]
        dataList = [float(x[0]) for x in dataList]
        gotScores = dict(zip(reqSeqs, dataList))

        # create final result from cache and webservice results
        scoreList = []
        for s in seqs:
            if s in gotScores:
                scoreList.append(gotScores[s])
            else:
                scoreList.append(scoreCache[s])

        scoreCache.update(gotScores)
        writeDict(scoreCache, cacheFname)
        return scoreList

    except urllib2.HTTPError, error:
        print("The request failed with status code: " + str(error.code))

        # Print the headers - they include the requert ID and the timestamp, which are useful for debugging the failure
        print(error.info())

        print(json.loads(error.read())) 

def trimSeqs(seqs, fiveFlank, threeFlank):
    """ given a list of 100bp sequences, return a list of sequences with the
    given number of basepairs 5' and 3' added from the middle position (pos 50) of
    the sequences
    """
    trimSeqs = []
    for s in seqs:
        trimSeqs.append(s[50+fiveFlank:50+threeFlank])
    return trimSeqs

def iterSvmRows(seqs):
    """ translate sequences to wang/sabatini/lander paper representation
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

def calcWangSvmScoresUsingR(seqs):
    """
    returns the SVM-calculated efficiency scores from the Wang/Sabatini/Lander paper
    This is using their original source code
    """
    writeSvmRows(seqs, "/tmp/temp.txt")
    wangSabDir = join(binDir, "src", "wangSabatiniSvm")
    cmd = "cd %s; R --slave --no-save -f scorer.R --args /tmp/temp.txt /tmp/temp.out" % wangSabDir
    print "running", cmd
    assert(os.system(cmd)==0)
    return parseSvmOut("/tmp/temp.out")


def calcAllScores(seqs, addOpt=[], doAll=False):
    """
    given 100bp sequences (50bp 5' of PAM, 50bp 3' of PAM) calculate all efficiency scores
    and return as a dict scoreName -> list of scores (same order).
    >>> calcAllScores(["CCACGTCTCCACACATCAGCACAACTACGCAGCGCCTCCCTCCACTCGGAAGGACTATCCTGCTGCCAAGAGGGTCAAGTTGGACAGTGTCAGAGTCCTG"], doAll=True)
    {'wang': [0.330643], 'chariRaw': [-0.15504833], 'chariRank': [54], 'doench': [0.03705739367422494], 'crisprScan': [0.391187032629], 'ssc': [-0.035894]}
    """
    scores = {}
    scores["wang"] = calcWangSvmScores(trimSeqs(seqs, -20, 0))
    scores["doench"] = calcDoenchScores(trimSeqs(seqs, -24, 6))
    scores["ssc"] = calcSscScores(trimSeqs(seqs, -20, 10))
    scores["crisprScan"] = calcCrisprScanScores(trimSeqs(seqs, -26, 9))

    chariScores = calcChariScores(trimSeqs(seqs, -20, 1))
    scores["chariRaw"] = chariScores[0]
    scores["chariRank"] = chariScores[1]

    if "fusi" in addOpt or doAll:
        scores["fusi"] = sendFusiRequest(trimSeqs(seqs, -24, 6)) # like doench
    return scores

def printScoreTabSep(seqs):
    " return tab-sep rows with all seqs "
    scoreDict = calcAllScores(seqs)
    scoreNames = scoreDict.keys()
    headers = ["fullSeq", "guideSeq"]
    headers.extend(scoreNames)

    print "\t".join(headers)
    for i, seq in enumerate(seqs):
        row = [seq, seq[30:53]]
        for scoreName in scoreNames:
            row.append(str(scoreDict[scoreName][i]))
        print "\t".join(row)

def test():
    global binDir
    binDir = "../crispor/bin"
    import doctest
    doctest.testmod()

if __name__=="__main__":
    if len(sys.argv)!=1:
        global binDir
        binDir = "../crispor/bin"
        inFname = sys.argv[1]
        seqs = open(inFname).read().splitlines()
        printScoreTabSep(seqs)
    else:
        test()
