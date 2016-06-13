# take all .ext.tab files in the effData directory, add the +/- 50bp context and write to 
# .context.tab files

# lots of checks to make sure the data is not screwed up in some way

from annotateOffs import *
import glob
import shutil

# very few guides have mapping problems
blackList = ['RPL36AL-4', 'Cd43-1']

def runBlatExtend(fname, seqs, db):
    """ extend 20mers to 100mers for eff score calculations, seqs is a dict id -> seq
    return a dict seqId -> list of (seq, genomePositionString)
    """
    # use cached BLAT results if possible
    # this is so BLAT can be run on a different server, with more RAM and the most 
    # sensitive settings 
    print fname
    blatBedFname = fname.replace(".guides.tab", ".blat.bed")
    print blatBedFname
    inNames = set()
    if isfile(blatBedFname):
        print "using %s instead of running BLAT" % blatBedFname
        shutil.copy(blatBedFname, "/tmp/temp.bed")
        for seqId, seq in seqs.iteritems():
            inNames.add(seqId)
    else:
        # write 19/20/23mers to fa file
        seqCount = 0
        ofh = open("/tmp/temp.fa", "w")
        for seqId, seq in seqs.iteritems():
            #if len(seq)==23:
                #seq = seq[:20]
            #if len(seq)==22:
                #seq = seq[:19]
            if len(seq) not in [20,22,23]:
                print ("Warning: sequence %s: %s is not 20 or 22 or 23 bp long. skipping" % (seqId, seq))
                continue
            ofh.write(">%s\n%s\n" % (seqId, seq))
            inNames.add(seqId)
            seqCount +=1 
        ofh.close()
        assert(seqCount > 0) # all sequences had wrong length

        print "running BLAT, writing to /tmp/temp.bed"
        if seqCount > 250 or db not in blatServers:
            #cmd = "blat -ooc=/gbdb/%s/11.ooc -dots=1 /gbdb/%s/%s.2bit /tmp/temp.fa /tmp/temp.psl -minScore=20 -noHead -stepSize=5 ; pslToBed /tmp/temp.psl /tmp/temp.bed" % (db, db, db)
            cmd = "blat -ooc=/gbdb/%s/11.ooc -dots=1000 /gbdb/%s/%s.2bit /tmp/temp.fa /tmp/temp.psl -minScore=20 -noHead -stepSize=5 ; pslToBed /tmp/temp.psl /tmp/temp.bed" % (db, db, db)
        else:
            blatServer, port = blatServers[db]
            cmd = "gfClient %s.soe.ucsc.edu %s /gbdb/%s /tmp/temp.fa /tmp/temp.psl -minScore=20 -nohead -minIdentity=100 -maxIntron=0 -dots=1 ; pslToBed /tmp/temp.psl /tmp/temp.bed" % (blatServer, port, db)
        print cmd
        assert(os.system(cmd)==0)

    matches = defaultdict(list) # seqId -> list of (chrom, start, end, strand)
    for line in open("/tmp/temp.bed"):
        chrom, start, end, name, score, strand = line.split()[:6]
        if db.startswith("hg") and ("_hap" in chrom or "random" in chrom or "chrUn" in chrom):
            continue
        start = int(start)
        end = int(end)
        # fix coordinates of 23bp long sequences
        if end-start==23:
            if strand=="+":
                end -= 3
            else:
                start += 3
        # fix coordinates of 22bp guide sequences (=19 bp guide)
        if end-start==22:
            if strand=="+":
                start -= 1
                end -= 3
            else:
                start += 3
                end += 1
        if (int(end)-int(start))!=20:
            continue
        matches[name].append( (chrom, int(start), int(end), strand) )

        # these matches cannot be found by gfClient so inject them here
        #if "Hamid_Guide_1a" in inNames:
            #matches["Hamid_Guide_1a"].append( ("chr3", 66933397, 66933420, "+") )
            #matches["Hamid_Guide_7b"].append( ("chr3", 66933396, 66933419, "-") )

    notFoundNames = inNames - set(matches)
    if len(notFoundNames)!=0:
        logging.warn("These sequences were not found with BLAT: %s" % ",".join(notFoundNames))

    filtMatches = {}
    for name, hits in matches.iteritems():
        if len(hits)!=1:
            logging.warn("dropping guide %s, matches multiple times: %s" % (name, hits))
        else:
            filtMatches[name] = hits

    assert(len(filtMatches)!=0) # no matches found ... is the genome database correct?

    # write extended matches to tempExt.bed file
    print "Creating /tmp/tempExt.bed with extended matches"
    ofh = open("/tmp/tempExt.bed", "w")
    positions = []
    for seqId, matchTuples in filtMatches.iteritems():
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
    print cmd
    os.system(cmd)
    seqs = parseFastaAsList(open("/tmp/tempExt.fa"))
    assert(len(seqs)==len(positions))

    ret = defaultdict(list)
    for seqData, pos in zip(seqs, positions):
        seqId, seq = seqData
        ret[seqId].append( (seq, pos) )

    return ret

def findAll(a_str, sub):
    start = 0
    while True:
        start = a_str.find(sub, start)
        if start == -1: return
        yield start
        #start += len(sub) # use start += 1 to find overlapping matches
        start += 1

def mapAndExtend20Mers():
    """
    code for new datasets:
    map 20mers to genome, extend to 100mer and write to context.tab file 
    """
    for fname in glob.glob("effData/*.guides.tab"):
        outFname = fname.replace(".guides.tab", ".context.tab")
        scoreFname = fname.replace(".guides.tab", ".scores.tab")
        if isfile(outFname):
            print "already there, %s" % outFname
            continue
        if isfile(scoreFname):
            print "already there, %s" % scoreFname
            continue

        print "working on %s" % fname
        # you can encode the db to map to in the filename
        dataset = basename(fname).split(".")[0]
        if "_" in fname:
            dataset, db = dataset.split("_")
        else:
            db = datasetToGenome[dataset]
        print "dataset %s: database for BLAT is %s" % (dataset, db)

        seqs = {}
        for row in iterTsvRows(fname):
            assert(row.guide not in seqs) # we don't want duplicated guides in here
            seqs[row.guide] = row.seq

        print "found %d sequences for BLAT" % len(seqs)
        matches = runBlatExtend(fname, seqs, db)

        notNggCount = 0
        fullRows = []
        for row in iterTsvRows(fname):
            for longSeq, pos in matches[row.guide]:
                longSeq = longSeq.upper()
                newRow = [row.guide, row.seq.upper(), db, pos, row.modFreq, longSeq]
                foundList =list(findAll(longSeq, newRow[1][:20]))
                if not 30 in foundList:
                    if 31 in foundList:
                        print("Warning: %s looks like a 19bp guide" % str(row))
                    else:
                        print "longSeq, newRow, foundList", longSeq, newRow, foundList
                        print (longSeq[:20].find(newRow[1]))
                        raise Exception("sequence extension did not work: %s" % str(newRow))

                if len(foundList)!=1:
                    print("Warning: extended sequence of %s contains the guide more than once" % newRow[0])

                if not (longSeq[51]=="G" and longSeq[52]=="G"):
                    print ("does not seem to have an NGG PAM, skipped: %s and PAM %s" % (str(newRow), longSeq[50:53]))
                    notNggCount +=1
                    continue

                fullRows.append(newRow)

        ofh = open(outFname, "w")
        writeRow(ofh, ["guide", "seq", "db", "pos", "modFreq", "longSeq"])
        for row in fullRows:
            writeRow(ofh, row)
        print "Wrote %s" % ofh.name
        print "Skipped %d guide because of missing NGG in the genome" % notNggCount
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

        # you can encode the db to map to in the filename
        if "_" in fname:
            db = basename(fname).split(".")[0].split("_")[1]
        else:
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

def extend20Mers():
    " mostly a copy of extend34Mers just for 20mers "
    for fname in glob.glob("effData/*.pos.tab"):
        outFname = fname.replace(".pos.tab", ".context.tab")
        if isfile(outFname):
            print "already there, %s" % outFname
            continue

        # db must be part of the filename e.g. dataset2016_hg19.pos.tab
        db = basename(fname).split(".")[0].split("_")[1]

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
            chrom, start, strand = row.chrom, row.start, row.strand
            start = int(start)
            end = int(start)+20
            seq = row.seq.upper()
            assert(len(seq)==20)

            if strand=="+":
                pamStart = end
            else:
                pamStart = start
            cxStart = pamStart - 50
            cxEnd = pamStart + 50
            newPos = "%s:%d-%d:%s" % (chrom, cxStart, cxEnd, strand)
            seqId = newPos
            bedRow = [chrom, str(cxStart), str(cxEnd), seqId, "0", strand]
            writeRow(bedFh, bedRow)
            # prep. final output row, except the long sequence which will be added later
            newRow = [row.guide, row.seq.upper(), db, newPos, row.modFreq]
            newRows.append(newRow)
            #if len(newRows)==100:
                #break

        # pull out the long sequences
        print "running twoBitToFa"
        bedFh.close()
        #cmd = "twoBitToFa /gbdb/%s/%s.2bit -bed=/tmp/tempExt.bed /tmp/tempExt.fa" % (db, db)
        cmd = "twoBitToFa /dev/shm/max/%s.2bit -bed=/tmp/tempExt.bed /tmp/tempExt.fa" % (db)
        os.system(cmd)
        extSeqs = parseFastaAsList(open("/tmp/tempExt.fa"))
        assert(len(extSeqs)==len(newRows))

        print "processing fa file"
        fullRows = []
        for (pos, extSeq), newRow in zip(extSeqs, newRows):
            #row, seq, db, chrom, start, strand, modFreq = newRow
            #print newRow, pos, extSeq, extSeq[50:70]
            seq = newRow[1]
            extSeq = extSeq.upper()
            newRow.append(extSeq)
            if not (extSeq.find(seq)==30):
                #print "extSeq", extSeq
                #print 'seq', seq, chrom, start, strand
                #print "find pos", extSeq.find(seq)
                print "sequence extension did not work: %s, extSeq %s, seq %s, pos %d" % (str(newRow), extSeq, seq, extSeq.find(seq))
                continue
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


# if we have an old .ext.tab file, then use it.
# this is data where I have manually select the best match for multi matching 20mers
# so we have the genome position already and I keep using it
extend34Mers()

# some new data comes as 20mers mapped to the genome: *.pos.tab files
# no need to run BLAT for these
extend20Mers()

# *.guides.tab
# new data comes as raw 20mers and have to be unique
# it's not mapped to the genome, so have to run BLAT
mapAndExtend20Mers()

