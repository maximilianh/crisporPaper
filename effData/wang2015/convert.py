from math import log

def getWangEssGenes():
    # get essential genes
    essGenes = set()
    for line in open("wang2015_SM_Table_S3.txt"):
        #Gene	sgRNAs included	KBM7 CS	KBM7 adjusted p-value	K562 CS	K562 adjusted p-value	Jiyoye CS	Jiyoye adjusted p-value	Raji CS	Raji adjusted p-value
        #A1BG	9	0.136	0.625949493	-0.220	0.951988369	0.096	0.630899802	-0.191	0.670542296
        #A1CF	7	-0.058	0.742293896	-0.984	0.100374265	0.217	0.386247018	-0.222	0.897305185
        if line.startswith("Gene"):
            continue
        fs = line.rstrip("\n").split("\t")
        gene = fs[0]
        p1, p2, p3, p4 = fs[3], fs[5], fs[7], fs[9]
        p1, p2, p3, p4 = float(p1), float(p2), float(p3), float(p4)
        if p1<0.01 and p2<0.01 and p3<0.01 and p4<0.01:
            essGenes.add(gene)
    return essGenes
    
def getHartEssGenes():
    # get essential genes
    essGenes = set()
    for line in open("../hart2016/essGenes/essential_sym_hgnc.csv"):
        fs = line.rstrip("\n").split("\t")
        gene = fs[0]
        essGenes.add(gene)
    return essGenes
    
# parse mod frequencies
modFreqs = {}
for line in open("aac7041_SM_Table_S2.txt"):
    if line.startswith("sgRNA"):
        continue
    fs = line.rstrip("\n").split("\t")
    name, kbm7Init1, kbm7Final1, kbmInit2, kbmFinal2 = fs[:5]
    vals = kbm7Init1, kbm7Final1, kbmInit2, kbmFinal2
    vals = [float(x) for x in vals]
    kbm7Init1, kbm7Final1, kbm7Init2, kbm7Final2 = vals
    # avg log fold change
    # lfc = log((kbm7Init1+kbmInit2+1)/2)-log((kbm7Final1+kbmFinal2+1)/2)
    #avgLfc = ((log(kbm7Final1+1)-log(kbm7Init1+1)) + (log(kbm7Final2+1)-log(kbm7Init2+1))) / 2.0
    avgLfc = (log(kbm7Final1+1, 2)-log(kbm7Init1+1, 2))
    #if name=="sgACTL6A_2":
        #print name, kbm7Init1, kbm7Final1, avgLfc
        #asdfdf
    modFreqs[name] = avgLfc

essGenes = getHartEssGenes()
    
ofh = open("../wang2015_hg19.pos.tab", "w")
ofh.write("guide\tseq\tchrom\tstart\tstrand\tmodFreq\n")
notFound = set()
for line in open("aac7041_SM_Table_S1.txt"):
    if line.startswith("sgRNA"):
        continue
    fs = line.rstrip("\n").split("\t")
    name, sym, chrom, start, strand, seq = fs[:6]
    if sym not in essGenes:
        continue
    if name not in modFreqs:
        notFound.add(name)
        continue
    modFreq = modFreqs[name]
    if start=="NA":
        continue
    start = str(int(start)-1)
    row = [name, seq, chrom, start, strand, str(modFreq)]
    ofh.write("\t".join(row))
    ofh.write("\n")

print "wrote to %s" % ofh.name
print "%d guides were not found in aac7041_SM_Table_S2.txt" % len(notFound)
