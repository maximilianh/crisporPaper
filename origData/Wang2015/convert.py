import glob
from collections import defaultdict

def compSeq(str1, str2):
    " return number of mismatches between str 1 and str 2"
    assert(len(str1)==len(str2))
    mm = 0
    for x, y in zip(str1, str2):
        if x!=y:
            mm += 1
    return mm

def revComp(seq):
    " rev-comp a dna sequence with UIPAC characters "
    revTbl = {'A' : 'T', 'C' : 'G', 'G' : 'C', 'T' : 'A', 'N' : 'N' , 'M' : 'K', 'K':'M'}
    newSeq = []
    for c in reversed(seq):
        newSeq.append(revTbl[c])
    return "".join(newSeq)

def parseSupp(fname, guideName):
    " return dict seq -> freq "
    firstSeq = None
    for line in open(fname):
        fs = line.strip().split()
        if len(fs)!=9:
            continue
        freq = fs[2]
        if freq=="UD":
            continue
        freq = float(freq)/100
        mmCount = fs[0]
        seq = fs[-1].upper()
        name = fs[1]
        if name.lower()=="mycn":
            # a typo in the suppl table sequence!!
            seq = "GAGGATGGGGAATGAGGAGTAGG"
        if firstSeq==None:
            firstSeq = seq
            mmCount = 0
            seqType = "on-target"
        else:
            mmCount = compSeq(firstSeq[:20], seq[:20])
            if mmCount>6:
                seq = revComp(seq)
                mmCount = compSeq(firstSeq[:20], seq[:20])
                assert(mmCount < 7)
            seqType = "off-target"

        #row = [guideName, seq.upper(), str(freq), seqType, name, str(mmCount)]
        row = [guideName, seq.upper(), str(freq), seqType]
        print "\t".join(row)

    
# get PCR + sequencing results frequencies from sup table 2 and 3
# the on-target frequencies were manually added here from supp fig 1
parseSupp("suppTable2.txt", "Wang_WAS-CR4")
parseSupp("suppTable3.txt", "Wang_WAS-CR5")
