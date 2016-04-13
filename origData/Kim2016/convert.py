# convert kim2016-suppTable3-hela.tsv to convert.tab

ofh = open("convert.tab", "w")
#ofh.write("name\tseq\tscore\ttype\n")

# has to look like this:
# name    seq     score   type
# Cho_c4bpb       AATGACCACTACATCCTCAAGGG 0.74678 on-target

# looks like this:
##name	Chromosome 	Location 	DNA seq at a Cleavage sites 	(-) RGEN 	(+) RGEN 	Validation
#VEGFA_OT	Chr6 	43737290 	GGGTGGGGGGAGTTTGCTCCAGG 	0.01% 	21.77% 	validated
#VEGFA1_02 	Chr15 	65637537 	GGaTGGaGGGAGTTTGCTCCTGG 	0.01% 	25.28% 	validated

digenomeSeqs = set(open("kimValidSeqs.txt").read().splitlines())
for s in digenomeSeqs:
    assert(len(s)==23)

for line in open("kim2016-suppTable3-hela.tsv"):
    if line.startswith("#"):
        continue
    fs = line.strip().split()
    fs = [s.strip() for s in fs]
    assert(fs[-1] in ["validated", "Invalidated"])
    if fs[-1] != "validated":
        continue
    #fs[0] = fs[0].replace("_", "-")
    if "OT" in fs[0]:
        seqType = "on-target"
    else:
        seqType = "off-target"
    fs[0] = fs[0].split("_")[0]
    val1 = str(float(fs[4].strip("%"))/100)
    val2 = str(float(fs[5].strip("%"))/100)
    if val2<val1:
        assert(False) # shouldn't happen, as they're "validated"
        #print "background too high"
        #continue
    seq = fs[3].upper()
    assert(len(seq)==23)
    if not seq in digenomeSeqs and seqType=="off-target":
        continue
    newRow = ["Kim16_%s" % fs[0], seq, val2, seqType]
    ofh.write("\t".join(newRow))
    ofh.write("\n")

print "wrote convert.tab"
