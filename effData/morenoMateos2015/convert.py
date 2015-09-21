ofh = open("../morenoMateos2015.ext.tab", "w")
ofh.write("guide\tseq\textSeq\tmodFreq\tposition\n")
for line in open("criprscan_input_moreno-mateos_vejnar.txt"):
    if line.startswith("#"):
        continue
    fs = line.strip("\n").split("\t")
    guide = fs[0]
    chrom, start, end = fs[1:4]
    seq = fs[4]
    strandStr = fs[5]

    strand = "+"
    if strandStr=="antisense":
        strand="-"

    #assert(strand==" or strand==-1)
    score = fs[-1]
    start, end = int(start), int(end)
    if strand=="+":
        start -= 4
        end += 7
    else:
        start -= 7
        end += 4
    posStr = "%s:%d-%d:%s" % (chrom, start, end, strand)
    row = [guide, seq+"NGG", "NNNN"+seq+"NGGNNNNNNN", score, posStr]
    ofh.write("\t".join(row))
    ofh.write("\n")
ofh.close()
print "Wrote to %s" % ofh.name
