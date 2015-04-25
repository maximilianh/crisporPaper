import glob

# had wang data as pseudo-bed files copied from Fig1, convert them here to 
# a simpler one-table format

def conv(guideName, mockFname, candFname, validFname):
    # get map name -> seq
    seqs = {}
    for line in open(mockFname):
        fs = line.rstrip("\n").strip().split()
        seq = fs[-1]
        name = fs[0]
        seqs[name] = seq

    # get map name -> frequency
    freqs = {}
    for line in open(validFname):
        fs = line.split()
        freqs[fs[0]] = fs[2]

    # go over digenome candidates
    for line in open(candFname):
        fs = line.rstrip("\n").strip().split()
        if line.lower().startswith("name"):
            continue
        name = fs[0]
        if name!="On-target":
            name = name.replace("-", "")

        pVal = fs[-1]
        if pVal=="N.A.":
            continue
        pVal = float(pVal)
        if pVal > 0.01:
            continue
        pVal = 1.0/pVal

        totalReads = fs[4]
        seqType = "off-target"
        if "On-target" in name:
            seqType = "on-target"
        if name=="name":
            continue
        #row = [guideName, seqs[name], str(pVal), freqs.get(name, "NA"), seqType]
        if name in freqs:
            freq = float(freqs[name].replace("%",""))/100
            row = [guideName, seqs[name], str(freq), seqType]
        print "\t".join(row)
        

#for fname in glob.glob("*.bed"):
    #for line in open(fname):
        ##fs = line.rstrip("\n").split()
        #chrom, start, end, seq, score = fs[:5]
        #fname = fname.split(".")[0]
        #seqType = "off-target"
        #if seq in targets:
            #seqType = "on-target"
        #row = [fname, seq, score, seqType]
        ##print "\t".join(row)

conv("Kim_"+"VEGF_A", "vegfaMock.txt", "vegfaCand.txt", "vegfaValid.txt")
conv("Kim_"+"HBB", "hbbMock.txt", "hbbCand.txt", "hbbValid.txt")
