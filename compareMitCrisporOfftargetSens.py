# list how many MIT is missing compared to crispor for each mismatch cutoff
from annotateOffs import *

# the two sequences ccr5-5 and ccr5-8 are the reverse complement of each other
# so they also containe two possible guide directions.
# this means that off-targets from MIT and crispor have artifically high mismatches, as the
# reverse strand will be considered an off-target. The MIT website unfortunately
# does not annotate the guide for an offtarget, so we don't know if the off-target
# relates to the forward or reverse strand guide.
# To avoid this problem, we simply skip these three guides for this analysis. 
skipGuides = ["Cho_ccr5-5", "Cho_ccr5-8", "Wang_WAS-CR5"]

def indexByMm(offDict, guideNames):
    """
    input: dict guideSeq -> list of off-target seqs
    output: dict mismatchCount -> list of off-target seqs
        and dict (mismatchCount, off-targetSeq) -> diffLogo
    given offtargets, return dict with mismatchCount -> list of offtarget seqs 
    REMOVES ALL OFFTARGET SEQUENCES THAT DON'T END with GG!!
    """
    mmToOffs = defaultdict(list)
    diffLogos = dict()
    for guideSeq, offs in offDict.iteritems():

        if guideSeq not in guideNames:
            print "skipping guide %s, not analyzed" % (guideSeq)
            continue
        guideName = guideNames[guideSeq]
        if guideName in skipGuides:
            print "skipping guide %s, to avoid revComp problems" % guideName
            continue

        for otSeq in offs:
            if not otSeq.endswith("GG"):
                continue
            mmCount, diffLogo = countMms(guideSeq[:20], otSeq[:20])
            if mmCount>=7:
                print guideName, repr(guideSeq), repr(otSeq), mmCount, diffLogo
                assert(False) # too many mismatches, sth is wrong
            diffLogos[(mmCount, otSeq)]=diffLogo
            mmToOffs[mmCount].append(otSeq)

    return mmToOffs, diffLogos

def parseCasOff(fname, guideSeqs):
    " parse casOff output file and return as nested dict guide -> otSeq -> None "
    shortToLong = {}
    for name, seq in guideSeqs.iteritems():
        shortToLong[seq[:20]] = seq

    print "parsing CasOff file %s" % fname

    ret = defaultdict(dict)
    notFoundGuides = set()
    for line in open(fname):
        if line.startswith("#"):
            continue
        fs = line.strip().split()
        guideSeq = fs[1][:20]
        if guideSeq not in shortToLong:
            notFoundGuides.add(guideSeq)
            continue
        guideSeq = shortToLong[guideSeq]
        otSeq = fs[2].upper()
        ret[guideSeq][otSeq] = 9999
    print "not found: guides %s" % (",".join(notFoundGuides))
    assert(len(ret)>0)
    return ret

def main():
    maxMismatches = 4
    guideValidOts, guideSeqs = parseOfftargets("out/annotFiltOfftargets.tsv", maxMismatches, False, ["GG", "AG", "GA"])
    casOffOffs = parseCasOff("casOffOfftargets/casOFFwithUPto4mm.txt", guideSeqs)
    mitOffs = parseMit("mitOfftargets", guideSeqs)
    crisporOffs = parseCrispor("crisporOfftargets", guideSeqs, 4)
    # inverse name -> seq to seq -> name mapping
    guideNames = {v: k for k, v in guideSeqs.items()}

    mm2Mit, _     = indexByMm(mitOffs, guideNames)
    mm2Crispor, crisprDiffLogos = indexByMm(crisporOffs, guideNames)
    mm2CasOff, _ = indexByMm(casOffOffs, guideNames)
    minMm = min(min(mm2Mit), min(mm2Crispor))
    maxMm = max(max(mm2Mit), max(mm2Crispor))

    ofh = open("out/mitCrisporSensDiff.tsv", "w")
    #row = [mm, len(mitOts), len(crispOts), len(casOffOts), len(crispNotMit), len(casOffNotMit), len(mitNotCrisp), len(casOffNotCrisp), mitMissSeq, crispMissSeqMit, crispMissSeqCasOff, ",".join(diffCounts)]
    headers = ["mismatches", "MIT_Predicted_NGG", "CRISPOR_Predicted_NGG", "CasOff_Predicted_NGG", "Crispor_notFoundMIT", "CasOff_notFoundMIT", "MIT_notFoundCrispor", "CasOff_notFoundCrispor", "MIT_notFoundCrispor_example", "Crispor_notFoundMIT_example", "CassOf_notFoundCrispor_example", "mismDistribution"]
    ofh.write("\t".join(headers)+"\n")

    for mm in range(minMm, maxMm+1):
        mitOts = set(mm2Mit.get(mm, []))
        crispOts = set(mm2Crispor.get(mm, []))
        casOffOts = set(mm2CasOff.get(mm, []))

        crispNotMit = crispOts - mitOts
        casOffNotMit = casOffOts - mitOts
        mitNotCrisp = mitOts - crispOts
        casOffNotCrisp = casOffOts - crispOts

        if len(crispNotMit)>0:
            mitMissSeq = list(crispNotMit)[0]
        else:
            mitMissSeq = "None"

        if len(mitNotCrisp)>0:
            crispMissSeqMit = list(mitNotCrisp)[0]
        else:
            crispMissSeqMit = "None"

        if len(casOffNotCrisp)>0:
            crispMissSeqCasOff = list(casOffNotCrisp)[0]
        else:
            crispMissSeqCasOff = "None"

        mmDiffLogos = [crisprDiffLogos[(mm, otSeq)] for otSeq in crispNotMit]
        
        # mismatches over length of missed off-targets, per position
        diffCounts = [0]*20
        for i in range(0,20):
            for logo in mmDiffLogos:
                if logo[i]=="*":
                    diffCounts[i]+=1
        diffCounts = [str(x) for x in diffCounts]

        row = [mm, len(mitOts), len(crispOts), len(casOffOts), len(crispNotMit), len(casOffNotMit), len(mitNotCrisp), len(casOffNotCrisp), mitMissSeq, crispMissSeqMit, crispMissSeqCasOff, ",".join(diffCounts)]
        row = [str(x) for x in row]
        ofh.write( "\t".join(row)+'\n')

    ofh.close()
    print "output written to %s" % ofh.name

main()
