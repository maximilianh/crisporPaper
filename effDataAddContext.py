# take all .ext.tab files in the effData directory, add the +/- 50bp context and write to 
# .context.tab files

# lots of checks to make sure the data is not screwed up in some way

from annotateOffs import *
import glob

# very few guides have mapping problems
blackList = ['RPL36AL-4', 'Cd43-1']

def runBlatExtend(seqs, db):
    """ extend 20mers to 100mers for eff score calculations, seqs is a dict id -> seq
    return a dict seqId -> list of (seq, genomePositionString)
    """
    # write 20mers to fa file
    inNames = set()
    ofh = open("/tmp/temp.fa", "w")
    for seqId, seq in seqs.iteritems():
        ofh.write(">%s\n%s\n" % (seqId, seq))
        inNames.add(seqId)
    ofh.close()

    print "running BLAT, writing to /tmp/temp.bed"
    blatServer, port = blatServers[db]
    #cmd = "gfClient %s.soe.ucsc.edu %s /gbdb/%s /tmp/temp.fa /tmp/temp.psl -minScore=20 -nohead -minIdentity=100 -maxIntron=0 -dots=1 ; pslToBed /tmp/temp.psl /tmp/temp.bed" % (blatServer, port, db)
    cmd = "blat /gbdb/%s/%s.2bit /tmp/temp.fa /tmp/temp.psl -minScore=20 -stepSize=4 -noHead  ; pslToBed /tmp/temp.psl /tmp/temp.bed" % (db, db)
    os.system(cmd)

    matches = defaultdict(list) # seqId -> list of (chrom, start, end, strand)
    for line in open("/tmp/temp.bed"):
        chrom, start, end, name, score, strand = line.split()[:6]
        if "_hap" in chrom or "random" in chrom or "chrUn" in chrom:
            continue
        if (int(end)-int(start))!=20:
            continue
        matches[name].append( (chrom, int(start), int(end), strand) )

    notFoundNames = inNames - set(matches)
    logging.warn("These sequences were not found with BLAT: %s" % ",".join(notFoundNames))

    for name, hits in matches.iteritems():
        if len(hits)!=1:
            logging.warn("Guide %s matches multiple times: %s" % (name, hits))

    # write extended matches to tempExt.bed file
    print "Creating /tmp/tempExt.bed with extended matches"
    ofh = open("/tmp/tempExt.bed", "w")
    positions = []
    for seqId, matchTuples in matches.iteritems():
        for matchTuple in matchTuples:
            chrom, start, end, strand = matchTuple
            if strand=="+":
                start = start - 30
                end = end + 50
            else:
                start = start - 50
                end = end + 30
            row = [chrom, str(start), str(end), seqId, "0", strand]
            ofh.write("\t".join(row)+"\n")
            positions.append( "%s:%d-%d:%s" % (chrom, start, end, strand))
    ofh.close()

    cmd = "twoBitToFa /gbdb/%s/%s.2bit -bed=/tmp/tempExt.bed /tmp/tempExt.fa" % (db, db)
    os.system(cmd)
    seqs = parseFastaAsList(open("/tmp/tempExt.fa"))
    assert(len(seqs)==len(positions))

    ret = defaultdict(list)
    for seqData, pos in zip(seqs, positions):
        seqId, seq = seqData
        ret[seqId].append( (seq, pos) )

    return ret

def mapAndExtend20Mers():
    """
    code for new datasets:
    map 20mers to genome, extend to 100mer and write to context.tab file 
    """
    for fname in glob.glob("effData/*.guides.tab"):
        outFname = fname.replace(".guides.tab", ".context.tab")
        if isfile(outFname):
            print "already there, %s" % outFname
            continue

        dataset = basename(fname).split(".")[0]
        db = datasetToGenome[dataset]

        seqs = {}
        for row in iterTsvRows(fname):
            assert(row.guide not in seqs)
            seqs[row.guide] = row.seq

        matches = runBlatExtend(seqs, db)

        fullRows = []
        for row in iterTsvRows(fname):
            for longSeq, pos in matches[row.guide]:
                newRow = [row.guide, row.seq.upper(), db, pos, row.modFreq, longSeq]
                if not (longSeq.find(newRow[1][:20])==30):
                    print (longSeq[:20].find(newRow[1]))
                    raise Exception("sequence extension did not work: %s" % str(newRow))
                if not (longSeq[51]=="G" and longSeq[52]=="G"):
                    raise Exception("not a NGG PAM? %s" % str(newRow))
                #newRow[1] = extSeq[30:53]
                fullRows.append(newRow)

        ofh = open(outFname, "w")
        writeRow(ofh, ["guide", "seq", "db", "pos", "modFreq", "longSeq"])
        for row in fullRows:
            writeRow(ofh, row)
        print "Wrote %s" % ofh.name
        #break

def extend34Mers():
    """
    extend legacy files with mapped 34mers
    """
    for fname in glob.glob("effData/*.ext.tab"):
        outFname = fname.replace(".ext.tab", ".context.tab")
        if isfile(outFname):
            print "already there, %s" % outFname
            continue

        db = datasetToGenome.get(basename(fname).split(".")[0], "hg19")
        print "Processing %s, db %s" % (fname, db)
        bedFh = open("/tmp/tempExt.bed", "w")
        newRows = []
        seenNames = set()
        for row in iterTsvRows(fname):
            if row.guide in seenNames:
                raise Exception("guide name %s already seen before" % row.guide)
            if row.guide in blackList:
                continue
            seenNames.add(row.guide)
            chrom, posStr, strand = row.position.split(":")
            start, end = posStr.split("-")
            start, end = int(start), int(end)
            seq = row.seq.upper()
            extSeq = row.extSeq.upper()

            assert(end-start)==34
            assert(len(seq)==23)
            assert(len(extSeq)==34)
            if not (extSeq.find(seq)==4):
                raise Exception("extended sequence does not include short seq at position 4: %s" % str(row))
            if strand=="+":
                pamStart = start+24
            else:
                pamStart = start + 10
            cxStart = pamStart - 50
            cxEnd = pamStart + 50
            newPos = "%s:%d-%d:%s" % (chrom, cxStart, cxEnd, strand)
            seqId = newPos
            bedRow = [chrom, str(cxStart), str(cxEnd), seqId, "0", strand]
            writeRow(bedFh, bedRow)
            # prep. final output row, except the long sequence which will be added later
            newRow = [row.guide, row.seq.upper(), db, newPos, row.modFreq]
            newRows.append(newRow)

        # pull out the long sequences
        bedFh.close()
        cmd = "twoBitToFa /gbdb/%s/%s.2bit -bed=/tmp/tempExt.bed /tmp/tempExt.fa" % (db, db)
        os.system(cmd)
        extSeqs = parseFastaAsList(open("/tmp/tempExt.fa"))
        assert(len(extSeqs)==len(newRows))

        fullRows = []
        for (pos, extSeq), newRow in zip(extSeqs, newRows):
            #print newRow, pos, extSeq, extSeq[50:70]
            extSeq = extSeq.upper()
            newRow.append(extSeq)
            if not (extSeq.find(newRow[1][:20])==30):
                print (extSeq[:20].find(newRow[1]))
                raise Exception("sequence extension did not work: %s" % str(newRow))
            if not (extSeq[51]=="G" and extSeq[52]=="G"):
                raise Exception("not a NGG PAM? %s" % str(newRow))
            newRow[1] = extSeq[30:53]
            fullRows.append(newRow)

        ofh = open(outFname, "w")
        writeRow(ofh, ["guide", "seq", "db", "pos", "modFreq", "longSeq"])
        for row in fullRows:
            writeRow(ofh, row)
        print "Wrote %s" % ofh.name
        #break

extend34Mers()
mapAndExtend20Mers()
