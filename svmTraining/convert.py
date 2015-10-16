ofh = open("chari.tab","w")
for line in open("chariSuppData1.txt"):
    if line.startswith("#"):
        continue
    fs = line.split("\t")
    if fs[0]!="":
        ofh.write("%s\t1\n" % fs[0])
    if fs[1]!="":
        ofh.write("%s\t-1\n" % fs[1])
ofh.close()
print "SVM training examples written to %s" % ofh.name


