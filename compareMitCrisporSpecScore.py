# compare guide specificity scores of MIT with crispor
from annotateOffs import *

def main():
    maxMismatches = 4
    guideValidOts, guideSeqs = parseOfftargets("annotFiltOfftargets.tsv", maxMismatches, False, ["GG", "AG", "GA"])
    mitOffs = parseMit("mitOfftargets", guideSeqs)
    crisporOffs = parseCrispor("crisporOfftargets", guideSeqs, 4)

    ofh = open("out/mitCrisporSpecDiff.tsv", "w")
    headers = ["guide", "MIT", "CRISPOR", "difference"]
    ofh.write("\t".join(headers)+"\n")
    for guideName, guideSeq in guideSeqs.iteritems():
        mitScore = calcMitGuideScore_offs(guideSeq, mitOffs[guideSeq])
        crisporScore = calcMitGuideScore_offs(guideSeq, crisporOffs[guideSeq])
        diff = mitScore - crisporScore
        row = [guideName, mitScore, crisporScore, diff]
        row = [str(x) for x in row]
        ofh.write( "\t".join(row)+'\n')
    ofh.close()

    print "output written to %s" % ofh.name

main()
