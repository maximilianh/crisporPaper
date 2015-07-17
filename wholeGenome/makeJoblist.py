# create a joblist file for parasol
import glob
ofh = open("jobList", "w")

for fname in glob.glob("jobFiles/*.fa"):
    fId = fname.split("/")[1].split(".")[0]
    cmd = "crispor.cgi hg19 {check in exists jobFiles/%s.fa} {check out exists out/guides/%s.tab} -o out/offs/%s.tab" % (fId, fId, fId)
    ofh.write(cmd+"\n")

print "wrote new file jobList"
