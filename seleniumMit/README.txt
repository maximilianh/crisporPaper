These scripts run 1000 randomly selected CDS refseq guide sequences through the 
mit CRISPR website and downloads the specificity scores.

randomLines ../wholeGenome/specScores.tab 1000 crisporSpecScores.random1000.refseqCDS.tab
python submitToMit.py crisporSpecScores.random1000.refseqCDS.tab 
# wait for a day
python getFromMit.py jobIds.txt
