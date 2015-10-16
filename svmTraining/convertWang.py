ofh = open("wang.tab","w")

# get dict name -> seq
# didn't work, no names found
#seqs = {}
#for line in open("wangS3.txt"):
    #if line.startswith("#"):
        #continue
    #fs = line.split("\t")
    #name = fs[0]
    #seq = fs[5]
    #seqs[name] = seq

for line in open("wangTraining.txt"):
    if line.startswith("#"):
        continue
    fs = line.split("\t")
    name = fs[0]
    vec = fs[6:-1]
    # convert vector to back to sequence
    assert(len(vec)==80)
    seq = []
    for i in range(0, 20):
        for x, nucl in enumerate("ACTG"):
            if vec[i*4+x]=="1":
                seq.append(nucl)
                break
    seq = "".join(seq)
    assert(len(seq)==20)

    score = float(fs[4])
    if score < -1.6: # from scorer.R received by Wang et al
        classStr = "1"
    else:
        classStr = "-1"

    ofh.write("%s\t%s\n" % (seq, classStr))
ofh.close()
print "Wang SVM training examples written to %s" % ofh.name


