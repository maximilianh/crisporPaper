import pickle, os, sys
sys.path.append("..")
from annotateOffs import parseFasta, revComp

scores = pickle.load(open("scores.pickle"))

strands = dict()
# create a new bed with extended coords
ofh = open("tempExt.bed", "w")
for line in open("temp.bed"):
    chrom, start, end, name, score, strand = line.split()[:6]
    # AAVS1-3 has matches on multiple chroms and even in the same locus
    if name=="AAVS1-3" and chrom!="chr19" and start!="55627615":
        continue
    if strand=="+":
       end = str(int(end)+3 )
    else:
       start = str(int(start)-3 )

    row = [chrom, start, end, name, score, strand]
    ofh.write("\t".join(row))
    ofh.write("\n")
    strands[name] = strand
ofh.close()

# get their seqs

headers = ["guide", "modFreq", "seq"]
print "\t".join(headers)

cmd = "twoBitToFa /gbdb/hg19/hg19.2bit -bed=tempExt.bed tempExt.fa"
os.system(cmd)
seqs = parseFasta(open("tempExt.fa"))

for seqId, seq in seqs.iteritems():
    score = scores[seqId]
    #if strands[seqId]=="-":
        #print seq
        #seq = revComp(seq)
        #print "after", seq
    row = [seqId, score, seq]
    print "\t".join(row)
