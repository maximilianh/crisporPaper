# filter the data from compEffScores for the plotHeat script

ofh = open("out/heatData.tsv", "w")
keepWords = ["Wang", "Doench MOLM", "Chari 293T", "Zebra", "Drosophila", "elegans", "Housden"]

doneLines = set()
for line in open("out/effScoreComp.tsv"):
    fs = line.split("\t")
    if line.startswith("dataset"):
        ofh.write(line)
        headers = fs
        continue
    name = line.split("\t")[0]

    if name in doneLines:
        continue
    doneLines.add(name)

    for word in keepWords:
        if word in name:
            ofh.write(line)
            break
    #total = sum([float(x) for x in fs[1:]])
    #print name, total
ofh.close()
print "output written to %s" % ofh.name
    

