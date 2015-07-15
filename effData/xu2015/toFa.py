import pickle
scores = dict()

for line in open("origData/Supplemental_Table_3.txt"):
    if line.startswith("#"):
        continue
    fs = line.rstrip("\n").split()
    name = "AAVS1"+"-"+fs[0]
    seq = fs[4]
    print ">%s\n%s" % (name, seq)
    score = fs[-1].replace("<2.0","2.0")
    scores[name] = score

for line in open("origData/Supplemental_Table_4.txt"):
    if line.startswith("#"):
        continue
    fs = line.rstrip("\n").split()
    score = fs[-1]
    name = fs[0]+"-"+fs[1]
    seq = fs[6]
    print ">%s\n%s" % (name, seq)
    score = fs[-1]
    scores[name] = score

pickle.dump(scores, open("scores.pickle", "w"))
