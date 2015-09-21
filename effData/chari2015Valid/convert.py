name2seq = {}
for line in open("chariS1.tab"):
    name, seq, casType = line.strip().split()[:3]
    if casType!="Sp":
        continue
    name2seq[name] = seq

lines = open("chari-S8.txt").read().splitlines()
cells = lines[0].strip().split()

ofhs = {}
for c in cells:
    ofh = open("../chari2015Valid_"+c.replace("-","")+".tab", "w")
    ofhs[c] = ofh
    ofh.write("guide\tseq\tmodFreq\n")

for l in lines[2:]:
    fs = l.strip().split()
    name = fs[0].split(".")[0]
    for i, cell in zip(range(3, len(fs), 3), cells):
        modFreq = float(fs[i])
        assert(modFreq < 100)
        seq = name2seq[name]
        ofhs[cell].write("%s\t%s\t%s\n" % (name, seq, modFreq))

print "wrote to %s" % ",".join([o.name for o in ofhs.values()])
