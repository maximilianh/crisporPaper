import gzip
i = 0
ofh = open("jobFiles/%05d.fa" % i, "w")

print "reading guides from allGuides.txt"
allGuides = set()
for line in gzip.open("allGuides.txt"):
    allGuides.add(line.strip())

print "writing to jobFiles/"
for guide in allGuides:
    if guide[:20].count("C") + guide[:20].count("G") >= 15:
        continue
    if guide[:20].count("A") + guide[:20].count("A") <= 4:
        continue

    i+=1
    if i%10==1:
        ofh.close()
        newFname = "jobFiles/%05d.fa" % (i/10)
        ofh = open(newFname, "w")
        #print "new file", newFname
    ofh.write(">%s\n%s\n" % (guide, guide))
ofh.close()
print i
