# run all guides from the input file through the MIT website

from selenium import webdriver
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.common.by import By
import time, sys

# read input file
inFname = sys.argv[1]
seqs = []
for line in open(inFname):
    seqs.append(line.split()[0])
seqs = set(seqs)

# don't resubmit anything that we already have
doneSeqs = set()
for line in open("jobIds.txt"):
    doneSeqs.add(line.split()[0])
todoSeqs = seqs - doneSeqs

print "%d seqs to do" % len(todoSeqs)

driver = webdriver.Firefox()
#driver = webdriver.Chrome()
driver.implicitly_wait(30) # seconds

email = "maximilianh@gmail.com"

ofh = open("jobIds.txt", "a")

for i, seq in enumerate(todoSeqs):
    if seq=="AAAGAGGAAGCTGTTTGTCAAGG":
        continue
    if seq=="GAGGCGGCCAGTGTAAGCGGAGG":
        continue
    print i, len(todoSeqs), seq
    driver.get("http://crispr.mit.edu")

    nameBox = driver.find_element_by_id("search-name")
    emailBox = driver.find_element_by_id("email-address")
    seqBox = driver.find_element_by_id("sequence_submission_area")

    nameBox.send_keys(seq)
    emailBox.send_keys(email)
    seqBox.send_keys(seq)
    nameBox.submit()

    driver.find_element_by_id("job-container")
    time.sleep(1)
    jobId = driver.current_url.split("/")[-1]
    ofh.write( seq+"\t"+jobId+"\n")
    ofh.flush()
    time.sleep(3)

driver.quit()

