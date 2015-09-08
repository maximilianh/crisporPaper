# take all .ext.tab files in the effData directory, add the +/- 50bp context and write to 
# .context.tab files

# lots of checks to make sure the data is not screwed up in some way

from annotateOffs import *
import glob

# very few guides have mapping problems
blackList = ['RPL36AL-4', 'Cd43-1']

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
        if not (extSeq.find(newRow[1])==30):
            print (extSeq.find(newRow[1]))
            raise Exception("sequence extension did not work: %s" % str(newRow))
        if not (extSeq[51]=="G" and extSeq[52]=="G"):
            raise Exception("not a NGG PAM? %s" % str(newRow))
        fullRows.append(newRow)

    ofh = open(outFname, "w")
    writeRow(ofh, ["guide", "seq", "db", "pos", "modFreq", "longSeq"])
    for row in fullRows:
        writeRow(ofh, row)
    print "Wrote %s" % ofh.name
    #break
        
