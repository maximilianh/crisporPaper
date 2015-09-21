targets = {}
for line in open("suppTable1.txt"):
    fs = line.strip().split()
    name = fs[1]
    id = fs[0]
    seq = fs[2]
    #targets[ name+"."+ id]= seq
    targets[seq] = 0.0

rows = []
guideSums = {}
countFh = open("readCounts.tab", "w")
countFh.write("otSeq\treadCount\n")

for line in open("SuppTable7.txt"):
    fs = line.strip().split()
    if len(fs)<5:
        continue
    if fs[0]=="Sample":
        continue
    #print fs
    guideName = fs[0].split("-")[0]

    seq = fs[1]+fs[2]
    freq = fs[5]
    readCount = fs[-3]
    guideSeq = seq[:20]
    if guideSeq in targets:
        sType = "on-target"
    else:
        sType = "off-target"

    countFh.write("%s\t%s\n" % (seq, readCount))

    #if float(freq)==0.0:
        #continue
    row = [ guideName, seq, freq, sType]
    rows.append(row)
    if not guideName in guideSums:
        guideSums[guideName] = 0.0
    guideSums[guideName]+= float(freq)
    #print targets

#print guideSums
for row in rows:
    #print "raw", float(row[2]), row
    # do not correct
    #row[2] = str( float(row[2]) / guideSums[row[0]] )
    row[0] = "Hsu_"+row[0]
    print "\t".join(row)

