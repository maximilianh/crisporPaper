# run all guides from specScores.rnd100.tab through the MIT website

from selenium import webdriver
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.common.by import By
import time, sys

# read input file
inFname = sys.argv[1]
jobIds = []
for line in open(inFname):
    jobIds.append(line.strip().split())

doneSeqs = set()
for line in open("seqScores.txt"):
    doneSeqs.add(line.split()[0])

driver = webdriver.Firefox()
driver.implicitly_wait(10) # seconds

ofh = open("seqScores.txt", "a")

for inSeq, jobId in jobIds:
    if inSeq[:20] in doneSeqs:
        continue
    driver.get("http://crispr.mit.edu/guides/%s" % jobId)

    qscoreEls = driver.find_elements_by_class_name("qscore")
    qscores = [q.text for q in qscoreEls]

    seqEls = driver.find_elements_by_class_name("guide")
    seqs = [s.text for s in seqEls]
    assert(len(seqs)>0)
    assert(len(qscores)>0)
    assert(len(qscores)==len(seqs))

    doneSeqs = set()
    for seq, qscore in zip(seqs, qscores):
        if seq[:20] != inSeq[:20]:
            continue
        if seq in doneSeqs:
            continue
        doneSeqs.add(seq)
        ofh.write( "%s\t%s\n" % (seq, qscore))
        ofh.flush()

driver.quit()

