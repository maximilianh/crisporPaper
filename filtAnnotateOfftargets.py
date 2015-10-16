# filter the big annotated offtarget list offtargets.tsv to out/offtargetsFilt.tsv
# - remove two outlier guides with a high GC content
# - remove all offtargets with a frequency < 0.001

# take this raw list of off-target sequences and add annotations to it, like
# number of mismatches, gc content, off-target score, efficiency scores, etc
# - output overview data on offtargets with likely "bulge" effects?
from annotateOffs import *

def makeOutRows(inRows, targetSeqs, specScores):
    rows = []
    for row in inRows:
        guideSeq = targetSeqs[row.name]
        otSeq = row.seq
        mmCount, diffLogo = countMms(guideSeq[:-3], otSeq[:-3])
        if len(guideSeq)==23:
            otScore = calcHitScore(guideSeq[:-3], otSeq[:-3])
        else:
            otScore = "NA_not20mer"
        guideGc = gcCont(guideSeq[:20])
        bulgeRnaMm, bulgeRnaGuideSeqs, bulgeRnaOtSeqs, bulgeRnaLogos = findGappedSeqs(guideSeq, otSeq, mmCount-3)
        bulgeDnaMm, bulgeDnaOtSeqs, bulgeDnaGuideSeqs, bulgeDnaLogos = findGappedSeqs(otSeq, guideSeq, mmCount-3)

        specScore = specScores[guideSeq]
        otRow = [row.name, guideSeq, specScore, otSeq, str(guideGc), \
                float(row.score), str(mmCount), otScore, diffLogo, \
                bulgeRnaMm, ",".join(bulgeRnaGuideSeqs), ",".join(bulgeRnaOtSeqs), \
                bulgeDnaMm, ",".join(bulgeDnaGuideSeqs), ",".join(bulgeDnaOtSeqs) \
                ]
        rows.append(otRow)

    rows.sort(key=operator.itemgetter(5), reverse=True)
    return rows

def annotateOfftargets(inFname, outFname, specScores):
    # take the raw list of off-target sequences and add annotations to it, like
    # number of mismatches, gc content, off-target score, efficiency scores, etc
    inRows, targetSeqs = parseRawOfftargets(inFname, removeCellLine=True)
    print "%d offtargets, %d guides" % (len(inRows), len(targetSeqs))

    rows = makeOutRows(inRows, targetSeqs, specScores)

    # write out rows
    headers = ["name", "guideSeq", "guideSpecScore4MM", "otSeq", "guideGc", "readFraction", "mismatches", "otScore", "diffLogo","bulgeRnaMmCount", "bulgeRnaGuideSeq", "bulgeRnaOtSeq", "bulgeDnaMmCount", "bulgeDnaGuideSeq", "bulgeDnaOtSeq"]
    ofh = open(outFname, "w")
    ofh.write( "\t".join(headers) )
    ofh.write( "\n")
    for row in rows:
        assert(len(row)==len(headers))
        row = [str(x) for x in row]
        ofh.write( "\t".join(row)+"\n")
    print "wrote %s, %d rows" % (ofh.name, len(rows))

def filterOfftargets(inFname, outFname):
    # filter the big annotated offtarget list offtargets.tsv to out/offtargetsFilt.tsv
    # - remove two outlier guides with a high GC content
    # - remove all offtargets with a frequency < 0.001
    # - output overview data on offtargets with likely "bulge" effects?

    gcConts = {}
    headers = None
    for row in iterTsvRows(inFname):
        if row.type == "on-target":
            rowGc = gcCont(row.seq[:20])
            gcConts[row.name] = rowGc
        headers = row._fields

    ofh = open(outFname, "w")
    ofh.write("\t".join(headers)+"\n")

    count = 0
    guides = set()
    guideSeqs = set()
    remGuides = set()
    #print gcConts
    filtNames =["Tsai_HEK293_sgRNA4", "Tsai_VEGFA_site2"]
    for row in iterTsvRows(inFname):
        if row.name in filtNames:
            continue
        #if gcConts[row.name]>=75:
            #remGuides.add(row.name)
            #continue
        if row.score < 0.001:
            continue
        ofh.write("\t".join(row)+"\n")
        count += 1
        guides.add(row.name)
        if row.type=="on-target":
            guideSeqs.add(row.seq)

    print "removed these guides: %s" % ",".join(filtNames)
    print "kept %d guide names (->tests of identical guides in diff. cells)" % len(guides)
    print "kept %d different guide sequences" % len(guideSeqs)
    print "kept %d off-targets" % count
    print "output written to %s" % ofh.name

def getSpecScores():
    " obtain spec scores by parsing off-target lists. Use a cache to speed up subsequent runs. "
    TMPFNAME = "/tmp/specScores.pickle"
    if isfile(TMPFNAME):
        print "reading spec scores  from temp file %s" % TMPFNAME
        return pickle.load(open(TMPFNAME))

    maxMismatches = 4
    crisporOffs = parseCrispor("crisporOfftargets", None, maxMismatches)
    specScores = {}
    for guideSeq in crisporOffs.keys():
        specScores[guideSeq] = calcMitGuideScore_offs(guideSeq, crisporOffs[guideSeq])
    pickle.dump(specScores, open(TMPFNAME, "w"))
    return specScores

def main():
    specScores = getSpecScores()

    inFname = "offtargets.tsv"
    annotFname = "out/annotOfftargets.tsv"
    print "annotating without filter:"
    annotateOfftargets(inFname, annotFname, specScores)

    filtFname = "out/offtargetsFilt.tsv"
    filtAnnotFname = "out/annotFiltOfftargets.tsv"

    filterOfftargets(inFname, filtFname)
    annotateOfftargets(filtFname, filtAnnotFname, specScores)

main()
