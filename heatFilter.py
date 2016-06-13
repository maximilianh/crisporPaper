# filter the data from compEffScores for the plotHeat script

import sys

# if True, do not remove housden and gc-based scores
doSuppl = False

ofh = open("out/heatData.tsv", "w")
if len(sys.argv)==2 and sys.argv[1]=="part2":
    # fig 5, only our datasets
    keepWords = ["Schoenig", "Eschstruth", "in vivo", "Shkumatava", "Concordet Hs/Mm/Rn"]
    sortOrder = ["Concordet Hs/Mm/Rn", "Eschstruth", "Schoenig", "In Vivo", "Shkumatava"]

elif len(sys.argv)==2 and sys.argv[1]=="part1":
    # fig 4, only published datasets, remove outliers
    keepWords = ["Wang/Xu HL60", "Ghandi", "Hart", "Koike", "Koike-Yusa/Xu", "Chari 293T", "Varshney",'Gagnon', 'Moreno-Mateos', "elegans", 'Hct1162', 'Doench 2016', 'Doench 2014 Mouse', 'Ren']
    # Doench 2014 Mouse
    sortOrder = ["Wang/Xu HL60", "Doench 2014", "Koike", "Chari", "Doench 2016", "Hart", "Ghandi", "Farboud", "Ren", "Varshney", "Gagnon", "Moreno-Mateos"]

else:
    # supp file, the whole enchilada
    keepWords = ["Wang/Xu HL60", "Ghandi", "Hart", "Koike", "Koike-Yusa/Xu", "Chari 293T", "Chari K562",  "Varshney",'Gagnon', 'Moreno-Mateos', "elegans", 'Hct1162', 'Doench 2014 MOLM13/NB4/TF1','Doench 2016', 'Doench 2014 Mouse', 'Ren', "Liu", "Housden", "Wang 2015", "Eschstruth", "Shkumatava Zebrafish", "In Vivo", "Schoenig K562 LacZ Rank", "Concordet Hs/Mm/Rn"]
    sortOrder = ["Wang/Xu HL60", 'Wang 2015', 'Doench 2014 MOLM13/NB4/TF1',"Doench 2014 Mouse", "Doench 2016", "Koike", "Chari 293T", "Chari K562", "Hart Rpe", "Hart Hct116-1 Lib 1", "Hart Hct116-2 Lib 1", "Ghandi", "Farboud", "Ren", "Liu", "Eschstruth", "Housden", "Varshney", "Gagnon", "Moreno-Mateos", "In Vivo"]
    doSuppl = True
    #keepWords = ["Liu", "Ghandi", "hart2016", "Wang", "wang", "Doench 2014 MOLM", "Doench 2014 Mouse", "Chari 293T", "Zebra", "Drosophila", "elegans", "Housden"]

maxCols = -2
if doSuppl:
    maxCols = None

doneLines = set()
outRows = []
for line in open("out/effScoreComp.tsv"):
    line = line.replace("Ghandi", "Gandhi")
    fs = line.split("\t")
    if line.startswith("dataset"):
        headers = fs
        remCol2 = headers.index("Housden Score")
        if not doSuppl:
            del fs[remCol2]
        ofh.write("\t".join(fs[:maxCols]))
        ofh.write("\n")
        continue
    name = line.split("\t")[0]

    if name in doneLines:
        continue
    doneLines.add(name)

    for word in keepWords:
        if word in name:
            #ofh.write(line)
            #del fs[remCol1]
            if not doSuppl:
                del fs[remCol2]
            #ofh.write("\t".join(fs[:maxCols]))
            #ofh.write("\n")
            outRows.append(fs[:maxCols])
            break
    #total = sum([float(x) for x in fs[1:]])
    #print name, total

# first take out rows in sort order
sortedRows = []
doneRows = []
for sortKey in sortOrder:
    for row in outRows:
        if sortKey in row[0]:
            sortedRows.append(row)
            break
# add the rest
for row in outRows:
    if row not in sortedRows:
        sortedRows.append(row)
# output
for row in sortedRows:
    ofh.write("\t".join(row))
    ofh.write("\n")
ofh.close()


print "output written to %s" % ofh.name
    

